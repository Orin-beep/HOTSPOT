# Test version (to optimize the model size)
import torch
import pickle as pkl
import torch.nn as nn


class SelfAttention(nn.Module):
    def __init__(self, embed_size, heads):
        super(SelfAttention, self).__init__()
        self.embed_size = embed_size
        self.heads = heads
        self.head_dim = embed_size // heads # 512/8=64

        assert (
            self.head_dim * heads == embed_size
        ), "Embedding size needs to be divisible by heads"

        self.values = nn.Linear(self.head_dim, self.head_dim, bias=False)   # No bias
        self.keys = nn.Linear(self.head_dim, self.head_dim, bias=False)
        self.queries = nn.Linear(self.head_dim, self.head_dim, bias=False)
        self.fc_out = nn.Linear(heads * self.head_dim, embed_size)  # 8*64 -> 512

    def forward(self, values, keys, query, mask):
        # Get number of training examples
        N = query.shape[0]

        value_len, key_len, query_len = values.shape[1], keys.shape[1], query.shape[1]  # 400
        # Split the embedding into self.heads different pieces
        values = values.reshape(N, value_len, self.heads, self.head_dim)    # (200, 400, 512)->(200,400,8,64)
        keys = keys.reshape(N, key_len, self.heads, self.head_dim)
        query = query.reshape(N, query_len, self.heads, self.head_dim)
        values = self.values(values)  # (N, value_len, heads, head_dim)
        keys = self.keys(keys)  # (N, key_len, heads, head_dim)
        queries = self.queries(query)  # (N, query_len, heads, heads_dim)
        # Einsum does matrix mult. for query*keys for each training example
        # with every other training example, don't be confused by einsum
        # it's just how I like doing matrix multiplication & bmm

        energy = torch.einsum("nqhd,nkhd->nhqk", [queries, keys])
        # queries shape: (N, query_len, heads, heads_dim),
        # keys shape: (N, key_len, heads, heads_dim)
        # energy: (N, heads, query_len, key_len)

        # Mask padded indices so their weights become 0
        if mask is not None:
            energy = energy.masked_fill(mask == 0, float("-1e20"))  # through softmax

        # Normalize energy values similarly to seq2seq + attention
        # so that they sum to 1. Also divide by scaling factor for
        # better stability
        attention = torch.softmax(energy / (self.embed_size ** (1 / 2)), dim=3)
        # attention shape: (N, heads, query_len, key_len)

        out = torch.einsum("nhql,nlhd->nqhd", [attention, values]).reshape(
            N, query_len, self.heads * self.head_dim
        )
        # attention shape: (N, heads, query_len, key_len)
        # values shape: (N, value_len, heads, heads_dim)
        # out after matrix multiply: (N, query_len, heads, head_dim), then
        # we reshape and flatten the last two dimensions.

        out = self.fc_out(out)
        # Linear layer doesn't modify the shape, final shape will be
        # (N, query_len, embed_size)

        return out


class TransformerBlock(nn.Module):
    def __init__(self, embed_size, heads, dropout, forward_expansion):
        super(TransformerBlock, self).__init__()
        self.attention = SelfAttention(embed_size, heads)
        self.norm1 = nn.LayerNorm(embed_size)
        self.norm2 = nn.LayerNorm(embed_size)
        self.feed_forward = nn.Sequential(
            nn.Linear(embed_size, forward_expansion * embed_size),  # hidden 4*512
            nn.ReLU(),
            nn.Linear(forward_expansion * embed_size, embed_size),
        )

        self.dropout = nn.Dropout(dropout)

    def forward(self, value, key, query, mask):
        attention = self.attention(value, key, query, mask)

        # Add skip connection, run through normalization and finally dropout
        x = self.dropout(self.norm1(attention + query))
        forward = self.feed_forward(x)
        out = self.dropout(self.norm2(forward + x)) # two residuals
        return out


class Encoder(nn.Module):
    def __init__(
        self,
        src_vocab_size,
        embed_size,
        num_layers,
        heads,
        device,
        forward_expansion,
        dropout,
        max_length,
    ):

        super(Encoder, self).__init__()
        self.embed_size = embed_size
        self.device = device
        self.word_embedding = nn.Embedding(src_vocab_size, embed_size)  # 512
        self.position_embedding = nn.Embedding(max_length, embed_size)

        self.layers = TransformerBlock(
                    embed_size,
                    heads,
                    dropout=dropout,
                    forward_expansion=forward_expansion,
                )

        self.dropout = nn.Dropout(dropout)

    def forward(self, x, mask):
        N, seq_length = x.shape
        positions = torch.arange(0, seq_length).expand(N, seq_length).to(self.device)
        out = self.dropout((self.word_embedding(x) + self.position_embedding(positions)))
        # In the Encoder the query, key, value are all the same, it's in the
        # decoder this will change. This might look a bit odd in this case.
        out = self.layers(out, out, out, mask)

        return out


class Transformer(nn.Module):
    def __init__(self,
        src_vocab_size,
        src_pad_idx,
        embed_size=512,
        num_layers=1,   # one encoder
        forward_expansion=4,
        heads=8,
        dropout=0,
        device="cpu",
        max_length=100,
        feat_size=64,
        output_num=2
    ):
        super(Transformer, self).__init__()

        self.encoder = Encoder(
            src_vocab_size,
            embed_size,
            num_layers,
            heads,
            device,
            forward_expansion,
            dropout,
            max_length,
        )

        self.src_pad_idx = src_pad_idx
        self.device = device
        self.fc1 = nn.Linear(max_length*embed_size, feat_size)
        self.bn = nn.BatchNorm1d(feat_size)

    def make_src_mask(self, src):
        src_mask = (src != self.src_pad_idx).unsqueeze(1).unsqueeze(2)  # 200, 1, 1, 400
        return src_mask.to(self.device)

    def forward(self, src):
        src_mask = self.make_src_mask(src)
        enc_src = self.encoder(src, src_mask)
        enc_src = enc_src.reshape(enc_src.shape[0], -1) # 200, 400*512
        feat = self.bn(self.fc1(enc_src))   # feature by transformer
        return feat


class hotspot_model(nn.Module): # hotspot: two transformers
    def __init__(self,
        src_vocab_size=2,
        embed_size=512,
        num_layers=1,   # one encoder
        forward_expansion=4,
        heads=8,
        dropout=0,
        device="cpu",
        output_num=2,
        dropout2=0
    ):
        super(hotspot_model, self).__init__()

        #self.feat_size1 = 160   # total 192 feature size
        #self.feat_size2 = 32
        self.feat_size1 = 106
        self.feat_size2 = 22
        self.src_pad_idx = 0
        self.max_length1 = 400  # MCL protein sentence length
        self.max_length2 = 50   # mob+mpf sentence length
        self.src_vocab_size1 = 108274+1 # num_pcs + zero padding
        self.src_vocab_size2 = 14   # 9 MOB types + 4 MPF types + 1 zero padding
        self.Trans1 = Transformer(self.src_vocab_size1,self.src_pad_idx,device=device,max_length=self.max_length1,dropout=dropout,feat_size=self.feat_size1, embed_size=256)    # MCL protein Transformer
        self.Trans2 = Transformer(self.src_vocab_size2,self.src_pad_idx,device=device,max_length=self.max_length2,dropout=dropout,feat_size=self.feat_size2, embed_size=16)     # MOB/MPF feature Transformer

        self.out = nn.Linear(self.feat_size1+self.feat_size2+22, output_num)
        self.dropout = nn.Dropout(dropout)


    def forward(self, src):
        sentence1 = src[:, :self.max_length1]
        sentence2 = src[:, self.max_length1:self.max_length1+self.max_length2]
        sentence3 = src[:, self.max_length1+self.max_length2:]
        feat1 = self.Trans1(sentence1)
        feat2 = self.Trans2(sentence2)
        feat=torch.hstack((feat1,feat2, sentence3))
        feat = self.dropout(feat)
        logit = self.out(feat)
        return logit, feat
