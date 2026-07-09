import torch 
import argparse
import os
import sys
import yaml 
from tqdm import tqdm
import json 
import numpy as np
import pandas as pd
sys.path.append(os.environ.get("SAFARI_PATH", "."))
from src.models.sequence.long_conv_lm import ConvLMHeadModel
from src.dataloaders.datasets.hg38_char_tokenizer import CharacterTokenizer

try:
    from tokenizers import Tokenizer  
except:
    pass

class HG38Encoder:
    "Encoder inference for HG38 sequences"
    def __init__(self, model_cfg, ckpt_path, max_seq_len, nlayer):
        self.max_seq_len = max_seq_len
        self.nlayer = nlayer
        self.model, self.tokenizer = self.load_model(model_cfg, ckpt_path, max_seq_len, nlayer)
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = self.model.to(self.device)

    def encode(self, seqs):
            
        results = []

        # sample code to loop thru each sample and tokenize first (char level)
        for seq in tqdm(seqs):
            
            if isinstance(self.tokenizer, Tokenizer):
                tokenized_seq = self.tokenizer.encode(seq).ids
            else:
                tokenized_seq = self.tokenizer.encode(seq)
            
            # can accept a batch, shape [B, seq_len, hidden_dim]
            logits, hidden_states = self.model(torch.tensor([tokenized_seq]).to(device=self.device))

            # Using head, so just have logits
            results.append(hidden_states)

        return results
        
            
    def load_model(self, model_cfg, ckpt_path, max_seq_len, nlayer):
        config = yaml.load(open(model_cfg, 'r'), Loader=yaml.FullLoader)
        config['model_config']['n_layer'] = nlayer 
        config['model_config']['layer']['l_max'] = max_seq_len + 2
        model = ConvLMHeadModel(**config['model_config'])
        
        state_dict = torch.load(ckpt_path, map_location='cpu')

        # loads model from ddp by removing prexix to single if necessary
        torch.nn.modules.utils.consume_prefix_in_state_dict_if_present(
            state_dict["state_dict"], "model."
        )

        model_state_dict = state_dict["state_dict"]

        # need to remove torchmetrics. to remove keys, need to convert to list first
        for key in list(model_state_dict.keys()):
            if "torchmetrics" in key:
                model_state_dict.pop(key)

        model.load_state_dict(state_dict["state_dict"])

        # setup tokenizer
        if config['tokenizer_name'] == 'char':
            print("**Using Char-level tokenizer**")

            # add to vocab
            tokenizer = CharacterTokenizer(
                characters=['A', 'C', 'G', 'T', 'N'],
                model_max_length=self.max_seq_len + 2,  # add 2 since default adds eos/eos tokens, crop later
		add_special_tokens=False,
            )
            print(tokenizer._vocab_str_to_int)
        else:
            raise NotImplementedError("You need to provide a custom tokenizer!")

        return model, tokenizer
        
        
if __name__ == "__main__":
    
    SAFARI_PATH = os.getenv('SAFARI_PATH', '.')

    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--model_cfg",
        default=f"{SAFARI_PATH}/configs/evals/hyena_small_150b.yaml",
    )
    
    parser.add_argument(
        "--ckpt_path",
        default=f"",
        help="Path to model state dict checkpoint"
    )

    parser.add_argument(
        "--seq_file",
        default=f"",
        help="Path to fasta file with sequences to encode"
    )
    
    parser.add_argument(
      "--output_file",
      default=f"",
      help="Path to output embeddings file"
    )
    
    parser.add_argument(
      "--max_seqlen",
      default=f"",
      help="Path to output embeddings file"
    )
    
    parser.add_argument(
      "--nlayers",
      default=f"",
      help="don't add two"
    )
    
    parser.add_argument(
      "--batch_size",
      default=f"",
      help="don't add two"
    )
        
    args = parser.parse_args()
        
    task = HG38Encoder(args.model_cfg, args.ckpt_path, max_seq_len=int(args.max_seqlen), nlayer=int(args.nlayers))
    print('Successfully loaded encoder.', flush =True)
    # sample sequence, can pass a list of seqs (themselves a list of chars)
    seqs = []
    names = []
    current_seq = []
    with open(args.seq_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    seqs.append(''.join(current_seq))
                    current_seq = []
                names.append(line[1:])
            else:
                current_seq.append(line)
    if current_seq:
        seqs.append(''.join(current_seq))
    # verify that the sequences are loaded correctly
    if not (len(seqs) == len(names)):
        print('Error: Number of sequences and names do not match.', flush=True)
        sys.exit(1)
    print('Successfully loaded data.', flush=True)
    
    # if seqs is too long for memory, can batch it
    # determine batch size based on length of all sequences
    total_seq_len = sum([len(seq) for seq in seqs])
    max_seq_len_per_batch = 400000
    batch_size =  int(args.batch_size)

    with open(args.output_file, "w") as f:
        for i in range(0, len(seqs), batch_size):
            print('Before encode.', flush=True)
            logits = task.encode(seqs[i:i+batch_size])
            print('After encode.',flush=True)
            #print(len(logits),flush=True)
            #print(logits[0].shape,flush=True)
            #print(logits[1].shape,flush=True)
            #print(len(logits[0]),flush=True)
            #print(logits[0].mean(dim=1).shape,flush=True)
            my_names = names[i:i+batch_size]
            print(f"{i+batch_size} sequences processed")
            for i in range(len(logits)):
                #print(logits[i].shape,flush=True)
                # get the mean of the entire hidden states
                hidden_states_mean = logits[i][0].mean(dim=0)
                out_list = [my_names[i]] + hidden_states_mean.tolist()
                f.write('\t'.join([str(x) for x in out_list]) + '\n')
    print('Successfully encoded sequences.', flush=True)
