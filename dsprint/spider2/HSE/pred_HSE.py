#!/usr/bin/env python
import numpy, os, sys

def read_pssm(pssm_file):
  # this function reads the pssm file given as input, and returns a LEN x 20 matrix of pssm values.

  # index of 'ACDE..' in 'ARNDCQEGHILKMFPSTWYV'(blast order)
  idx_res = (0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18)
  
  # open the two files, read in their data and then close them
  fp = open(pssm_file, 'r')
  lines = fp.readlines()
  fp.close()

  # declare the empty dictionary with each of the entries
  aa = []
  pssm = []
  
  # iterate over the pssm file and get the needed information out
  for line in lines:
    split_line = line.split()
    # valid lines should have 32 points of data.
    # any line starting with a # is ignored
    if (len(split_line) == 44) and (split_line[0] != '#'):
      aa_temp = split_line[1]
      aa.append(aa_temp)
      pssm_temp = [-float(i) for i in split_line[2:22]]
      pssm.append([pssm_temp[k] for k in idx_res])
  
  return aa, pssm

def read_hhm(hhm_file):
  f = open(hhm_file)
  line=f.readline()
  while line[0]!='#':
      line=f.readline()
  f.readline()
  f.readline()
  f.readline()
  f.readline()
  seq = []
  extras = numpy.zeros([0,10])
  prob = numpy.zeros([0,20])
  line = f.readline()
  while line[0:2]!='//':
      lineinfo = line.split()
      seq.append(lineinfo[0])  
      probs_ = [2**(-float(lineinfo[i])/1000) if lineinfo[i]!='*' else 0. for i in range(2,22)]
      prob = numpy.concatenate((prob,numpy.matrix(probs_)),axis=0)
      
      line = f.readline()
      lineinfo = line.split()
      extras_ = [2**(-float(lineinfo[i])/1000) if lineinfo[i]!='*' else 0. for i in range(0,10)]
      extras = numpy.concatenate((extras,numpy.matrix(extras_)),axis=0)
      
      line = f.readline()
      assert len(line.strip())==0
      
      line = f.readline()
  #return (''.join(seq),prob,extras)
  return (seq,numpy.concatenate((prob,extras),axis=1))

def read_spider2(spider2_file):
  # This function reads the spider 2 file given as input. 

  # open spider 2 file and read data.
  fp = open(spider2_file,'r')
  lines = fp.readlines()
  fp.close()
  
  # declare the empty dictionary with each of the entries
  aa = []
  ss = []
  asa = []
  ttpp = []

  # iterate over the file and get the needed information out
  dict_ASA0 = dict(zip("ACDEFGHIKLMNPQRSTVWY", (115, 135, 150, 190, 210, 75, 195, 175, 200, 170,185, 160, 145, 180, 225, 115, 140, 155, 255, 230)))
  b1 = not lines[0].startswith('#index')
  for k,line in enumerate(lines):
    split_line = line.split()
    if (len(split_line) != 11) or (split_line[0][0] == '#'): continue
    # valid lines should have 11 points of data.
    # any line starting with a # is ignored
    if b1:
         aa_temp = split_line[1]
         aa.append(aa_temp)
         ss_temp = [float(split_line[i]) for i in (-3,-2,-1)]
         ss.append(ss_temp)
         asa.append([float(split_line[3])/dict_ASA0[aa_temp]])
                 
         tt_temp = numpy.array([numpy.radians(float(split_line[i])) for i in (6,7)])
         tt_temp = numpy.concatenate((numpy.sin(tt_temp), numpy.cos(tt_temp)), axis=0)/2 + 0.5
         pp_temp = numpy.array([numpy.radians(float(split_line[i])) for i in (4,5)])
         pp_temp = numpy.concatenate((numpy.sin(pp_temp), numpy.cos(pp_temp)), axis=0)/2 + 0.5
         ttpp.append(numpy.concatenate((tt_temp, pp_temp), axis=0))

    else:
         aa_temp = split_line[1]
         aa.append(aa_temp)
         ss_temp = [float(i) for i in split_line[3:6]]
         ss.append(ss_temp)
         asa.append([float(split_line[6])])
         ttpp_temp = [float(i) for i in split_line[7:11]]
         tt_temp = numpy.array([numpy.radians(float(i)) for i in split_line[7:9]])
         tt_temp = numpy.concatenate((numpy.sin(tt_temp), numpy.cos(tt_temp)), axis=1)/2 + 0.5
         pp_temp = numpy.array([numpy.radians(float(i)) for i in split_line[9:11]])
         pp_temp = numpy.concatenate((numpy.sin(pp_temp), numpy.cos(pp_temp)), axis=1)/2 + 0.5
         ttpp.append(numpy.concatenate((tt_temp, pp_temp), axis=1))

  return aa, ss, asa, ttpp

def get_input(pssm_file, spider2_file, hhm_file='NULL'):
  # this function takes a path to a pssm file and finds the pssm + phys 7 input
  # to the NN in the required order - with the required window size (8).
  
  
  # define the dictionary with the phys properties for each AA
  phys_dic = {'A': [-0.350, -0.680, -0.677, -0.171, -0.170, 0.900, -0.476],
              'C': [-0.140, -0.329, -0.359, 0.508, -0.114, -0.652, 0.476],
              'D': [-0.213, -0.417, -0.281, -0.767, -0.900, -0.155, -0.635],
              'E': [-0.230, -0.241, -0.058, -0.696, -0.868, 0.900, -0.582],
              'F': [ 0.363, 0.373, 0.412, 0.646, -0.272, 0.155, 0.318],
              'G': [-0.900, -0.900, -0.900, -0.342, -0.179, -0.900, -0.900],
              'H': [ 0.384, 0.110, 0.138, -0.271, 0.195, -0.031, -0.106],
              'I': [ 0.900, -0.066, -0.009, 0.652, -0.186, 0.155, 0.688],
              'K': [-0.088, 0.066, 0.163, -0.889, 0.727, 0.279, -0.265],
              'L': [ 0.213, -0.066, -0.009, 0.596, -0.186, 0.714, -0.053],
              'M': [ 0.110, 0.066, 0.087, 0.337, -0.262, 0.652, -0.001],
              'N': [-0.213, -0.329, -0.243, -0.674, -0.075, -0.403, -0.529],
              'P': [ 0.247, -0.900, -0.294, 0.055, -0.010, -0.900, 0.106],
              'Q': [-0.230, -0.110, -0.020, -0.464, -0.276, 0.528, -0.371],
              'R': [ 0.105, 0.373, 0.466, -0.900, 0.900, 0.528, -0.371],
              'S': [-0.337, -0.637, -0.544, -0.364, -0.265, -0.466, -0.212],
              'T': [ 0.402, -0.417, -0.321, -0.199, -0.288, -0.403, 0.212],
              'V': [ 0.677, -0.285, -0.232, 0.331, -0.191, -0.031, 0.900],
              'W': [ 0.479, 0.900, 0.900, 0.900, -0.209, 0.279, 0.529],
              'Y': [ 0.363, 0.417, 0.541, 0.188, -0.274, -0.155, 0.476],
              'X': [ 0.0771,-0.1536, -0.0620, -0.0762, -0.1451,  0.0497, -0.0398],
              'Z': [ 0.0771,-0.1536, -0.0620, -0.0762, -0.1451,  0.0497, -0.0398]}
  
  pssm_aa, pssm = read_pssm(pssm_file)
  spider2_aa, ss, asa, ttpp = read_spider2(spider2_file)
  
  if pssm_aa != spider2_aa:
    print('PSSM residue sequence - ')
    print(pssm_aa)
    print('SPIDER2 residue sequence - ')
    print(spider2_aa)
    print("ERROR: residue sequence in pssm input and spider2 input are not the same.")
    sys.exit() 
  
  if hhm_file=='NULL':
    hhm = 0
    hhm_aa = 0
  else:
    hhm_aa, hhm = read_hhm(hhm_file)
    if pssm_aa != hhm_aa:
      print('PSSM residue sequence - ')
      print(pssm_aa)
      print('HHM residue sequence - ')
      print(hhm_aa)
      print("ERROR: residue sequence in hhm input is not the same as that in pssm and spider2 inputs.")
      sys.exit() 
     
  # calculate the length
  aa = pssm_aa
  length = len(aa)
    
  # set the phys7 data.
  phys = [phys_dic[i] for i in aa]
   
  return(pssm, aa, phys, ss, asa, ttpp, length, hhm)  

def window(feat, winsize=8):
  # apply the windowing to the input feature
  feat = numpy.array(feat)
  output = numpy.concatenate([numpy.vstack([feat[0]]*winsize), feat])
  output = numpy.concatenate([output, numpy.vstack([feat[-1]]*winsize)])
  output = [numpy.ndarray.flatten(output[i:i+2*winsize+1]).T for i in range(0,feat.shape[0])]
  return output
  
def window_data(*feature_types):
  n = len(feature_types[0])
  features = numpy.empty([n,0])

  for feature_type in feature_types:
    test = numpy.array(window(feature_type))
    features = numpy.concatenate((features, test), axis=1)
  
  return features
  

def sigmoid(input): 
  # apply the sigmoid function
  output = 1 / (1 + numpy.exp(-input))
  return(output)

def nn_feedforward(nn, input):
  input = numpy.matrix(input)

  # find the number of layers in the NN so that we know how much to iterate over
  num_layers = nn['n'][0][0][0][0]
  # num_input is the number of input AAs, not the dimentionality of the features
  num_input = input.shape[0]
  x = input

  # for each layer up to the final
  for i in range(1,num_layers-1):
    # get the bais and weights out of the nn
    W = nn['W'][0][0][0][i-1].T
    temp_size = x.shape[0]
    b = numpy.ones((temp_size,1))
    x = numpy.concatenate((b, x),axis=1)
    # find the output of this layer (the input to the next)
    x = sigmoid(x * W)

  # for the final layer.
  # note that this layer is done serpately, this is so that if the output nonlinearity
  # is not sigmoid, it can be calculated seperately.
  W = nn['W'][0][0][0][-1].T
  b = numpy.ones((x.shape[0],1))
  x = numpy.concatenate((b, x),axis=1)
  output = x*W
  pred = sigmoid(output)

  return pred

def load_NN(nn_filename):
  # load in the NN mat file.
  import scipy.io
  mat = scipy.io.loadmat(nn_filename)
    
  nn = mat['nn']
   
  return nn

  
def main(pssm_file, spider2_file, hhm_file, dat_nn1, dat_nn2, out_file):

  SS_order = ('C' 'E' 'H')

  pssm, aa, phys, ss, asa, ttpp, length, hhm = get_input(pssm_file, spider2_file, hhm_file)
  
  ## DO FIRST PREDICTIONS
  if type(hhm)==int:  # HHM is int when we don't pass in HHM file input. (hhm=0)
    input_feature = numpy.concatenate((numpy.ones((length,1))*length, window_data(pssm, phys, ss, asa, ttpp)), axis=1)
  else:
    input_feature = numpy.concatenate((window_data(hhm), numpy.ones((length,1))*length,window_data(pssm, phys, ss, asa, ttpp)), axis=1)
    #input_feature = window_data(hhm), length, window_data(pssm, phys, ss, asa, ttpp)
  norm_max = dat_nn1['high'][0][0][0]
  norm_min = dat_nn1['low'][0][0][0]
  # normalise the features
  hse_input_feature = (input_feature - numpy.tile(norm_min, (input_feature.shape[0], 1))) / numpy.tile((norm_max - norm_min), (input_feature.shape[0],1))
  # run NN feedforward
  pred_hse_1 = nn_feedforward(dat_nn1, hse_input_feature)
  pred_hse_1[:,0] *= 85
  pred_hse_1[:,1] *= 50
  pred_hse_1[:,2] *= 65
   
  # print first predictions to a file.
#  fp = open(pssm_file+'.hse1', 'w')
#  fp.write('#index\tAA\tCN\tHSEd\tHSEu\n')
#  for ind in range(0,len(aa)):
#    fp.write('%i\t%c\t%.3f\t%.3f\t%.3f\n' % (ind+1, aa[ind], pred_hse_1[ind,0], pred_hse_1[ind,1], pred_hse_1[ind,2]))
#  fp.close()
   
  ## DO SECOND PREDICTIONS  
  if type(hhm)==int:  # HHM is int when we don't pass in HHM file input. (hhm=0)
    input_feature = numpy.concatenate((numpy.ones((length,1))*length, window_data(pssm, phys, ss, asa, ttpp, pred_hse_1)), axis=1)
    #input_feature = numpy.concatenate((window_data(pred_hse_1), input_feature), axis = 1)
  else:
    input_feature = numpy.concatenate((window_data(hhm), numpy.ones((length,1))*length,window_data(pssm, phys, ss, asa, ttpp, pred_hse_1)), axis=1)
    #input_feature = numpy.concatenate((window_data(pred_hse_1), input_feature), axis = 1)
#
  norm_max = dat_nn2['high'][0][0][0]
  norm_min = dat_nn2['low'][0][0][0]
  # normalise the features
  hse_input_feature = (input_feature - numpy.tile(norm_min, (input_feature.shape[0], 1))) / numpy.tile((norm_max - norm_min), (input_feature.shape[0],1))
  # run NN feedforward
  pred_hse_2 = nn_feedforward(dat_nn2, hse_input_feature)
  pred_hse_2[:,0] *= 85
  pred_hse_2[:,1] *= 50
  pred_hse_2[:,2] *= 65
 
  # print predictions to a file.
  fp = open(out_file, 'w')
  fp.write('#index\tAA\tCN\tHSEu\tHSEd\n')
  for ind in range(0,len(aa)):
    fp.write('%i\t%c\t%.1f\t%.1f\t%.1f\n' % (ind+1, aa[ind], pred_hse_2[ind,0], pred_hse_2[ind,1], pred_hse_2[ind,2]))
  fp.close()

  return

  # DO ASA PREDICTION
  input_feature = window_data(pred_hse_2, pssm, phys, ss, asa, ttpp)
  nn = load_NN(nndir+'ASA_1.mat')
  norm_max = nn['high'][0][0][0]
  norm_min = nn['low'][0][0][0]
  # normalise the features
  hse_input_feature = (input_feature - numpy.tile(norm_min, (input_feature.shape[0], 1))) / numpy.tile((norm_max - norm_min), (input_feature.shape[0],1))
  # run NN feedforward
  pred_asa = nn_feedforward(nn, hse_input_feature)

  # print predictions to a file
  fp = open(pssm_file+'.asa1', 'w')
  fp.write('#index\tAA\tASA\n')
  for ind in range(0,len(aa)):
    fp.write('%i\t%c\t%.3f\n' % (ind+1, aa[ind], pred_asa[ind]))
  fp.close()

if __name__ == "__main__":
  # if there is no filename for the features to be written to given as input, don't write it
  if len(sys.argv) < 4:
    print >>sys.stderr, "Usage: RUN nndir pssm_dir hhm_dir spider2file... [-hsb]\n"
    sys.exit()
   
  nndir, pssm_dir, hhm_dir = sys.argv[1:4]
  dat_nn1 = load_NN(os.path.join(nndir, 'HSE_1.mat'))
  dat_nn2 = load_NN(os.path.join(nndir, 'HSE_2.mat'))

  suffix = 'hsa'
  if '-hsb' in sys.argv: suffix = 'hsb'
  if hhm_dir != 'NULL': suffix = 'h' + suffix

  fhhm = 'NULL'
  for x in sys.argv[4:]:
    if x[0] == '-': continue
    bn = os.path.basename(x)
    if bn.endswith('.spd3'): bn = bn[:-5]
    fpssm = os.path.join(pssm_dir, bn+'.pssm')
    if hhm_dir != 'NULL': fhhm = os.path.join(hhm_dir, bn+'.hhm')
    out_file = '%s.%s2' % (bn, suffix)
    main(fpssm, x, fhhm, dat_nn1, dat_nn2, out_file)
