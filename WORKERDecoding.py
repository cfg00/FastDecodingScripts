def run_single_fov_decoding(save_folder,fov,set_='',lib_fl = r'C:\Users\carlos\NMERFISH\codebook_AJBB_51hybe_DNA_V2_blank.csv'):
    dec = decoder_simple(save_folder,fov,set_=set_)
    if not os.path.exists(dec.decoded_fl):
        dec.get_XH(dec.fov,set_,ncols=3,nbits=51,th_h=5000,filter_tag = '')#number of colors match 
        dec.XH = dec.XH[dec.XH[:,-4]>0.25] ### keep the spots that are correlated with the expected PSF for 60X
        dec.load_library(lib_fl,nblanks=-1)
        dec.ncols = 3
        get_intersV2(dec,nmin_bits=3,dinstance_th=2,enforce_color=True,enforce_set=17*3,redo=False)
        get_icodesV3(dec,nmin_bits=3,iH=-3) ### saves a decoded....npz
        print("Decoding completed for {}.".format(dec.fov))

import sys
sys.path.append(r'C:\Users\carlos\NMERFISH')
from ioMicro import *
#from ioMicro import *
all_flds = glob.glob(r'\\merfish12.ucsd.edu\Public\Quan\D111_12_01_2023__AJBBqzbb_chr14\H*_DMER_*')
all_flds += glob.glob(r'\\merfish13.ucsd.edu\Public\Quan\D111_12_01_2023__AJBBqzbb_chr14\H*_DMER_*')
all_flds = np.array(all_flds)[np.argsort([get_iH(fld) for fld in all_flds])]
print("Found folders:",len(all_flds))
fov = 'Conv_zscan__200.zarr'
save_folder = r'\\merfish12.ucsd.edu\Public\Quan\D111_12_01_2023__AJBBqzbb_chr14\AnalysisDeconvolveCG_fixed'
dec = decoder_simple(save_folder,fov,set_='')
def load_segmentation(dec,segmentation_folder = r'\\merfish13.ucsd.edu\Public\Quan\D111_12_01_2023__AJBBqzbb_chr14\cell_segmentation',segm_tag='H1_DMER_1'):
    fl_segm = segmentation_folder+os.sep+dec.fov.replace('.zarr','')+'--'+segm_tag+'--CYTO_segm.npz'
    segm,shape = np.load(fl_segm)['segm'],np.load(fl_segm)['shape']
    #segm_ = resize(segm,shape)
    dec.im_segm_ = segm
    dec.shape=shape
    dec.segm_tag=segm_tag
scores_ref_fl = save_folder+os.sep+r'scoresRef.pkl'
def main_f_fov(save_folder =save_folder,fov='Conv_zscan__050',set_ = '',ncols=3,
          scores_ref_fl=scores_ref_fl,th=-5,force=False):
    save_fld_cell = os.path.dirname(save_folder)+os.sep+'best_per_cellComb'
    if not os.path.exists(save_fld_cell): os.makedirs(save_fld_cell)
    save_fl = save_fld_cell+os.sep+fov+'__XHfs_finedrft.npz'
    if not os.path.exists(save_fl) or force:
        dec = decoder_simple(save_folder,fov,set_)
        dec.save_fl=save_fl
        dec.ncols = ncols
        dec.load_decoded()

        #apply_fine_drift(dec,plt_val=True)
        #dec.save_folder= r'C:\Users\cfg001\Desktop\WTC11\flat_field'
        #apply_flat_field(dec,tag='Scope4_med_col_raw')
        #dec.save_folder= save_folder
        dec.ncols = ncols
        if False:
            #### prefiltering
            Hlog_min = np.log(np.nanmin(dec.XH_pruned[...,-3],axis=1))
            keep = Hlog_min>0
            cor_min = np.nanmin(dec.XH_pruned[...,-4],axis=1)
            keep&=cor_min>0.5
            cms = np.nanmean(dec.XH_pruned[...,:3],axis=1)[:,np.newaxis]
            dcenters = np.linalg.norm(cms-dec.XH_pruned[...,:3],axis=-1)
            dcenters_max = np.nanmax(dcenters,axis=-1)
            keep&=dcenters_max<1.5
            dec.XH_pruned,dec.icodesN=dec.XH_pruned[keep],dec.icodesN[keep]
        
        
        #scoresRefT = get_score_per_color(dec)
        scoresRefT = pickle.load(open(scores_ref_fl,'rb'))
        dec.dist_best = np.load(dec.decoded_fl)['dist_best']
        get_score_withRef(dec,scoresRefT,plt_val=True,gene=None,iSs = None,th_min=-7.5,include_dbits=True)
        dec.th=th
        plot_statistics(dec)
    
        #threshold the combined EM score
        keep = dec.scoreA>dec.th
        dec.XH_prunedf,dec.icodesNf=dec.XH_pruned[keep],dec.icodesN[keep]
        nbits = dec.XH_prunedf.shape[1]
        dec.XH_prunedF = np.concatenate([dec.XH_prunedf,np.repeat(dec.icodesNf,nbits).reshape(-1,nbits)[:,:,np.newaxis]],axis=-1)
        ### get the drift - to correct to the segmentation space
        dec.dic_drift = get_dic_drift(dec)
        
        load_segmentation(dec)
        dec.im_segm_ = expand_segmentation(dec.im_segm_, nexpand=5)
        
        ### Compute drift between the segmentation file and the reference drift file
        tzxy_seg = dec.dic_drift[dec.segm_tag] #np.round([dic_drift[key][0] for key in dic_drift if cp.segm_tag in key]).astype(int)
        ### Augment the fitting data with cell id
        resc = dec.im_segm_.shape/dec.shape
        XH_ = dec.XH_prunedF.copy()
        XH_[:,:,:3] = XH_[:,:,:3]-tzxy_seg ### bring fits to cell segmentation space - modified to -
        XC = (np.nanmean(XH_[:,:,:3],axis=1)*resc).astype(int) #rescale to segmentation size
        dec.XC = XC
        keep = np.all(XC>=0,axis=-1)&np.all(XC<dec.im_segm_.shape,axis=-1)
        icells = np.zeros(len(XC))
        icells[keep]=dec.im_segm_[tuple(XC[keep].T)]
        nbits = XH_.shape[1]
        icells = np.repeat(icells,nbits).reshape(-1,nbits)[:,:,np.newaxis]
        XH_f = np.concatenate([XH_,icells],axis=-1)
        dec.XH_f=XH_f
        XH_fs=keep_best_per_cell_fast(XH_f,nbest=20)
        
        np.savez(dec.save_fl,XH_fs=XH_fs)
        return dec
def get_dic_drift(dec):
    drifts,flds,fov_,fl_ref = np.load(dec.drift_fl,allow_pickle=True)
    return {os.path.basename(fld):drft[0] for drft,fld in zip(drifts,flds)}
def main_f(fov,try_mode=True):
    set_=''
    
    
    if True:#not os.path.exists(fl):
        if try_mode:
            try:
                compute_drift_V2(save_folder,fov,all_flds,set_='',redo=False,gpu=True)
                run_single_fov_decoding(save_folder,fov,set_='',lib_fl = r'C:\Users\carlos\NMERFISH\codebook_AJBB_51hybe_DNA_V2_blank.csv')
                dec = main_f_fov(fov=fov,th=-7,force=False)
            except:
                print("Failed within the main analysis:",fov)
        else:
            compute_drift_V2(save_folder,fov,all_flds,set_='',redo=False,gpu=True)
            run_single_fov_decoding(save_folder,fov,set_='',lib_fl = r'C:\Users\carlos\NMERFISH\codebook_AJBB_51hybe_DNA_V2_blank.csv')
            dec = main_f_fov(fov=fov,th=-7,force=False)
    
    
    return fov
    

from multiprocessing import Pool, TimeoutError    
if __name__ == '__main__':
    # start 4 worker processes
    fovs = [os.path.basename(fl) for fl in np.sort(glob.glob(all_flds[0]+os.sep+'*.zarr'))]
    item = fovs[56]
    main_f(item,try_mode=False)
    if True:
        with Pool(processes=7) as pool:
            print('starting pool')
            result = pool.map(main_f, fovs)
#conda activate cellpose2&&python C:\Users\carlos\NMERFISH\WorkerDecodingD111.py