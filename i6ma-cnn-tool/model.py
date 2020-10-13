
import numpy as np
import pandas as pd
import gc
import os
from itertools import *
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import to_categorical
from werkzeug.utils import secure_filename
from Bio import SeqIO

### Loading Physicochemical Properties
di_prop = pd.read_csv('di_selected.csv')
di_nucleotide = di_prop.iloc[:, 0].tolist()
pp_di = {}

for i in range(len(di_nucleotide)):
	pp_di[di_nucleotide[i]] = di_prop.iloc[i, 1:].tolist()
 
### Numerical Representations
lookup_table = []
for p in product('ATGC', repeat=2):
	w = ''.join(p)
	lookup_table.append(w)

### Monomer Binary Encoding
def monomer_matrix(seqs):
	X = np.empty((len(seqs), len(seqs[0])))
	alphabet = 'ATGC'

	for i in range(len(seqs)):
		for j in range(len(seqs[0])):
			X[i, j] = next((k for k, letter in enumerate(alphabet) if letter == seqs[i][j]))

	X = to_categorical(X, num_classes=4)

	return X
### Dimer Binary Encoding
def dimer_index_mapping(seqs):
	X = np.empty((len(seqs), len(seqs[0])-1))

	for i in range(len(seqs)):
		for j in range(len(seqs[0])-1):
			w = seqs[i][j:j+2]
			X[i,j] = lookup_table.index(w)

	X = to_categorical(X, num_classes=16)
	return X

### Dinucleotide Structural Properties
def di_prop(seqs):
	X = np.empty([len(seqs), len(seqs[0])-1, 9], dtype=float)

	for i in range(len(seqs)):
		for j in range(len(seqs[0])-1):
			word = seqs[i][j:j+2]
			value = pp_di[word]
			X[i, j] = value
	return X

### PseAAC Inspired
def pseaac(seqs):
	X = np.empty((len(seqs), len(seqs[0]), 4))

	ring_structure = { 'A':1, 'C':0, 'T':0, 'G':1 }
	functional_group = { 'A':1, 'C':1, 'T':0, 'G':0 }
	hydrogen_bonding = { 'A':0, 'C':1, 'T':0, 'G':1 }

	for i in range(len(seqs)):
		count_A = 0
		count_C = 0
		count_T = 0
		count_G = 0
		for j in range(len(seqs[0])):
			if seqs[i][j] == 'A':
				count_A += 1
				X[i, j, 3] = count_A / float(j+1)
			elif seqs[i][j] == 'C':
				count_C += 1
				X[i, j, 3] = count_C / float(j+1)
			elif seqs[i][j] == 'T':
				count_T += 1
				X[i, j, 3] = count_T / float(j+1)
			elif seqs[i][j] == 'G':
				count_G += 1
				X[i, j, 3] = count_G / float(j+1)

			X[i, j, 0] = ring_structure[seqs[i][j]]
			X[i, j, 1] = functional_group[seqs[i][j]]
			X[i, j, 2] = hydrogen_bonding[seqs[i][j]]

	return X
def change_file_name(file_name):
    new_name = 'static/uploads/dataset.txt'
    os.rename(file_name, new_name)
    return new_name

def handle_fa(filename):
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    seq = []
    names = []
    i = 0
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if len(sequence) != 41:
            assert False, 'Each sequence must have a length of 41nt.\nSequence {} has length {}nt'.format(i, len(sequence))
        seq.append(sequence)
        names.append(name)
        i+=1
    return seq,names

def prediction(file_name):
    '''
    fasta_file = change_file_name(file_name)
    # Handling FASTA files
    dataset = get_data(open(fasta_file, encoding='unicode_escape'), desc=True)
    seqs = []
    names = []

    for i in range(len(dataset)):
        if len(dataset[i].seq) != 41:
            assert False, 'Each sequence must have a length of 41nt.\nSequence {} has length {}nt'.format(i, len(dataset[i].seq))
        seqs.append(dataset[i].seq)
        names.append(dataset[i].name)
    '''
    seqs,names = handle_fa(file_name)
    ### Loading Pre-trained Model
    model = load_model('N6_final_model.h5')

    ### seqs -> is a list that contains dna sequences (upper case letters)
    batch_size = 1024
    pred = []
    n_rounds = np.ceil(len(seqs) / batch_size).astype('int')

    for i in range(n_rounds):
        batch_seqs = seqs[i*batch_size:i*batch_size+batch_size]

        X1 = monomer_matrix(batch_seqs)
        X2 = dimer_index_mapping(batch_seqs)
        X3 = di_prop(batch_seqs)
        X4 = pseaac(batch_seqs)

        pred.append(np.where(model.predict( [X1, X2, X3, X4], batch_size=batch_size) >= 0.5, 'Positive', 'Negative').ravel() )
        
        del X1, X2, X3, X4, batch_seqs
        gc.collect()

    pred_df = pd.DataFrame(columns=['Identifier','Sequence','Prediction'])
    pred = list(chain.from_iterable(pred))
    pred_df['Identifier'] = names
    pred_df['Sequence'] = seqs
    pred_df['Prediction'] = pred
    pred_df.to_csv('prediction.csv', index=False)
