import os.path
import glob
import dsprint.spider2.misc.pred_pssm as pred_pssm
import dsprint.spider2.HSE.pred_HSE as pred_HSE


try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <pssm_folder> <output_folder>')
        sys.exit(0)

    PSSM_FOLDER, OUTPUT_FOLDER = sys.argv[1:]
else:
    PSSM_FOLDER = snakemake.input.pssm_folder
    OUTPUT_FOLDER = snakemake.output.output_folder

if __name__ == '__main__':

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    nndir = os.path.join(os.path.dirname(pred_pssm.__file__), '..', 'dat/')
    dict1_nn = pred_pssm.load_NN(nndir + 'pp1.npz')
    dict2_nn = pred_pssm.load_NN(nndir + 'pp2.npz')
    dict3_nn = pred_pssm.load_NN(nndir + 'pp3.npz')
    list_nn = (dict1_nn, dict2_nn, dict3_nn)
    list_params = [list_nn, False]  # False = don't generate files at each iteration, just the last one

    for pssm_file in glob.glob(f'{PSSM_FOLDER}/*.pssm'):
        basename = os.path.splitext(os.path.basename(pssm_file))[0]
        aa, pssm = pred_pssm.read_pssm(pssm_file)
        spd_file = os.path.join(OUTPUT_FOLDER,  f'{basename}.spd')
        pred_pssm.pred1(list_params, aa, pssm, spd_file)

        spd3_file = os.path.join(OUTPUT_FOLDER, f'{basename}.spd3')
        for hsx in ('hsa', 'hsb'):
            hsx_folder = os.path.join(os.path.dirname(pred_HSE.__file__), f'{hsx}_full')
            dat_nn1 = pred_HSE.load_NN(os.path.join(hsx_folder, 'HSE_1.mat'))
            dat_nn2 = pred_HSE.load_NN(os.path.join(hsx_folder, 'HSE_2.mat'))
            hsx_file = os.path.join(OUTPUT_FOLDER, f'{basename}.{hsx}2')

            pred_HSE.main(pssm_file, spd3_file, 'NULL', dat_nn1, dat_nn2, hsx_file)
