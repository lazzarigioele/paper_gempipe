
import json 
from pathlib import Path
import glob
import pandas as pnd

from Bio import SeqIO, SeqRecord, SeqFeature, Seq

import gempipe
import cobra
cobra_config = cobra.Configuration()
cobra_config.solver = "cplex"

import multiprocess  # import multiprocessing
multiprocess.set_start_method('fork')  # starts a fresh python intepreter



def fasta_to_modeled(model, infiles, outfile):
    
    modeled_gids = set([g.id for g in model.genes]) -set(['spontaneous'])
    modeled_gids = set([gid[2:] if gid.startswith('G_') else gid for gid in modeled_gids])
    
    sr_list = []
    added = set()
    
    for infile in infiles: 
        for record in SeqIO.parse(infile, "fasta"):
            gid = '-'
            for attrib in record.description.split(' '):
                if attrib.startswith('[locus_tag='):
                    gid = attrib.replace('[locus_tag=', '').replace(']', '')
            if gid != '-' and gid in modeled_gids:
                sr = SeqRecord.SeqRecord(record.seq, id=gid, description='')
                sr_list.append(sr)
                added.add(gid)
    with open(outfile, 'w') as w_handler:
        count = SeqIO.write(sr_list, w_handler, "fasta")
        
    
    print(f'Recovered {len(added)} / {len(modeled_gids)}.')
    print("Missing", len(modeled_gids-added), ':', modeled_gids-added)
    return added



def gb_to_modeled(model, infiles, outfile):
        
    modeled_gids = set([g.id for g in model.genes]) -set(['spontaneous'])
    modeled_gids = set([gid[2:] if gid.startswith('G_') else gid for gid in modeled_gids])
        
    sr_list = []
    added = set()

    for infile in infiles: 
        # use SeqIO.read() if you expect only one record in the GenBank file.
        # alternatively, use SeqIO.parse() if there are multiple records in the file.
        records = SeqIO.parse(infile, 'genbank')
        for record in records:

            sequence_id = record.id
            sequence_description = record.description
            sequence_seqlen = len(record.seq)

            # access features:
            for feature in record.features:
                if feature.type == 'CDS' or feature.type == 'tRNA':
                    try: location = feature.location
                    except: location = None
                    try: gene_symbol = feature.qualifiers['gene'][0]
                    except: gene_symbol = None
                    try: aa_seq = feature.qualifiers['translation'][0]
                    except: aa_seq = None
                    if aa_seq == None:
                        aa_seq = str(feature.extract(record.seq).translate())  # reverse complement 
                    try: nt_seq = str(feature.extract(record.seq))
                    except: nt_seq = None
                    try: locus_tag = feature.qualifiers['locus_tag'][0]
                    except: locus_tag = None
                    try: old_locus_tag1 = feature.qualifiers['old_locus_tag'][0]
                    except: old_locus_tag1 = None
                    try: old_locus_tag2 = feature.qualifiers['old_locus_tag'][1]
                    except: old_locus_tag2 = None
                    try: protein_id = feature.qualifiers['protein_id'][0]
                    except: protein_id = None
                    
                    
                    for locus_tag in [locus_tag, old_locus_tag1, old_locus_tag2]:
                        if locus_tag in added:
                            continue
                        if locus_tag in modeled_gids:
                            if outfile.endswith('faa') and aa_seq != None:
                                sr = SeqRecord.SeqRecord(Seq.Seq(aa_seq), id=locus_tag, description='')
                                sr_list.append(sr)
                            if outfile.endswith('fna') and nt_seq != None:
                                sr = SeqRecord.SeqRecord(Seq.Seq(nt_seq), id=locus_tag, description='')
                                sr_list.append(sr)
                            added.add(locus_tag)

    with open(outfile, 'w') as w_handler:
        count = SeqIO.write(sr_list, w_handler, "fasta")

        
    print(f'Recovered {len(added)} / {len(modeled_gids)}.')
    print("Missing", len(modeled_gids-added), ':', modeled_gids-added)
    return added



def apply_medium(model):
    
    # reset all the exchanges:
    for r in model.reactions: 
        if len(r.metabolites)==1 and list(r.metabolites)[0].id.rsplit('_', 1)[-1]=='e':
            r.bounds = (0, 1000)

    # main inorganic
    model.reactions.get_by_id('EX_h_e').lower_bound = -1000
    model.reactions.get_by_id('EX_h2o_e').lower_bound = -1000
    model.reactions.get_by_id('EX_o2_e').lower_bound = -1000
    model.reactions.get_by_id('EX_co2_e').lower_bound = -1000
    
    # main sources:
    model.reactions.get_by_id('EX_glc__D_e').lower_bound = -10   
    model.reactions.get_by_id('EX_pi_e').lower_bound = -1000       
    model.reactions.get_by_id('EX_so4_e').lower_bound = -1000      
    model.reactions.get_by_id('EX_nh4_e').lower_bound = -1000      
    
    # trace element
    #model.reactions.get_by_id('EX_cbl1_e').lower_bound = -1000
    model.reactions.get_by_id('EX_ca2_e').lower_bound = -1000
    model.reactions.get_by_id('EX_cl_e').lower_bound = -1000
    model.reactions.get_by_id('EX_cobalt2_e').lower_bound = -1000
    model.reactions.get_by_id('EX_cu2_e').lower_bound = -1000
    model.reactions.get_by_id('EX_fe2_e').lower_bound = -1000
    model.reactions.get_by_id('EX_fe3_e').lower_bound = -1000
    #model.reactions.get_by_id('EX_tungs_e').lower_bound = -1000
    model.reactions.get_by_id('EX_k_e').lower_bound = -1000
    model.reactions.get_by_id('EX_mg2_e').lower_bound = -1000
    model.reactions.get_by_id('EX_mn2_e').lower_bound = -1000
    model.reactions.get_by_id('EX_zn2_e').lower_bound = -1000
    model.reactions.get_by_id('EX_mobd_e').lower_bound = -1000
    #model.reactions.get_by_id('EX_na1_e').lower_bound = -1000
    #model.reactions.get_by_id('EX_ni2_e').lower_bound = -1000
    
    
    
def save_recipe(model, flavour='gempipe', medium_id='medium', outpath='./'):
    # recipe for gapseq are defined by hand     
    
    rows = [] if flavour=='carveme' else {'name': medium_id, 'exchanges': {}}
    
    # iterate exchange reactions:
    for r in model.reactions: 
        if len(r.metabolites)==1 and list(r.metabolites)[0].id.rsplit('_', 1)[-1] == 'e':
            if r.bounds[0] < 0: 
                
                puremid = list(r.metabolites)[0].id   # remove 'EX_' and '_e'
                if puremid.startswith('EX_'): puremid = puremid[3:]
                if puremid.endswith('_e'): puremid = puremid[:-2]
                
                if flavour=='carveme': 
                    rows.append({'medium': medium_id, 'description': '-', 'compound': puremid, 'name': '-'})
                else: 
                    if flavour=='bactabolize' and (r.id in ['EX_co2_e', 'EX_o2_e']) :
                        continue  # the --atmosphere_type aerobic option will be used from the command line
                    rows['exchanges'][r.id] = r.lower_bound
    
    if flavour=='carveme': 
        rows = pnd.DataFrame.from_records(rows)
        rows.to_csv(outpath + medium_id + '.tsv' , sep='\t', index=False)
    else :
        with open(outpath + medium_id + '.json', "w") as json_file:
            json.dump(rows, json_file)
            
            
            
def convert_fasta_to_genbank_with_cds(fasta_file, outfolder):

    fasta_sequences = SeqIO.parse(fasta_file, "fasta")
    genbank_records = []

    # convert each FASTA sequence to a GenBank record with a single CDS features:
    for index, seq_record in enumerate(fasta_sequences, start=1):

        genbank_record = SeqRecord.SeqRecord(
            Seq.Seq(str(seq_record.seq)),
            id=seq_record.id,
            name=seq_record.id,
            description=seq_record.id,
            annotations={"molecule_type": "protein"}
        )

        cds_feature = SeqFeature.SeqFeature(
            SeqFeature.FeatureLocation(start=0, end=len(seq_record)),
            type="CDS",
            qualifiers={
                "translation": str(seq_record.seq),
                "product": seq_record.id,
                "protein_id": seq_record.id,
                "locus_tag": seq_record.id
            }
        )

        genbank_record.features.append(cds_feature)
        genbank_records.append(genbank_record)

    SeqIO.write(genbank_records, outfolder + '/' + Path(fasta_file).stem + '.gb', "genbank")
    
    
    
def gapfilling_for_bactabolize(ref_file, file, outfolder, dataset='01_klebsiella'):
    
    ref_model = cobra.io.read_sbml_model(ref_file)
    model = cobra.io.load_json_model(file)
    
    apply_medium(ref_model)
    if dataset=='02_ralstonia':
        ref_model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0   
        ref_model.reactions.get_by_id('EX_gln__L_e').lower_bound = -10
    if dataset=='03_pseudomonas':
        ref_model.reactions.get_by_id('EX_na1_e').lower_bound = -1000   
        ref_model.reactions.get_by_id('EX_ni2_e').lower_bound = -1000
        
    ref_res = ref_model.optimize()
    ref_obj = ref_res.objective_value
    ref_status = ref_res.status
    
    apply_medium(model)
    if dataset=='02_ralstonia':
        model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0   
        model.reactions.get_by_id('EX_gln__L_e').lower_bound = -10
    if dataset=='03_pseudomonas':
        model.reactions.get_by_id('EX_na1_e').lower_bound = -1000   
        model.reactions.get_by_id('EX_ni2_e').lower_bound = -1000  

    res = model.optimize()
    obj = res.objective_value
    status = res.status

    print('Gapfilling', file, 'ref_obj', ref_obj, 'ref_status', ref_status, 'obj', obj, 'status', status)  # log
    
    if obj < 0.6 or status == 'infeasible':

        rids = gempipe.perform_gapfilling(model, ref_model, minflux= 0.6, nsol=1, verbose=True)
        print('Gapfilling', file, 'with', rids, gempipe.get_solver(ref_model), gempipe.get_solver(model))  # log
        if rids != None:   # at least a solution was found
            for rid in rids: 
                gempipe.import_from_universe(model, ref_model, rid)

    cobra.io.save_json_model(model, outfolder + '/' + Path(file).stem + '.json')
        
        
        
def simulate_biolog(file, outfolder, seed=False, dataset='01_klebsiella', tool='carveme'): 
    
    if file.endswith('.json'):
        model = cobra.io.load_json_model(file)
    if file.endswith('.xml'):
        model = cobra.io.read_sbml_model(file)
    
    if not seed:   # medium is already applied in gapseq models
        apply_medium(model)
        
        if dataset == '03_pseudomonas':
            if tool in ['gempipe', 'bactabolize']:   # satisfy iJN1463 requirements in ref-based recons
                model.reactions.get_by_id('EX_na1_e').lower_bound = -1000   
                model.reactions.get_by_id('EX_ni2_e').lower_bound = -1000 
    
    print('Biolog', file, 'with', model.slim_optimize(), gempipe.get_solver(model))  # log
    biolog_preview = gempipe.biolog_preview(model, seed=seed)
    
    biolog_preview.to_csv(outfolder + '/' + Path(file).stem + '.csv')
    
    
    
def get_biolog_C_mappings():
    
    # import the gempipe's biolog mappings
    biolog_mappings = gempipe.get_biolog_mappings()

    # from the biolog_mappings, keep only PM1 and PM2: 
    interesting_subs = []
    for substrate, row in biolog_mappings.iterrows():
        for well in eval(row['PM']):
            if well.startswith('PM01:') or well.startswith('PM02:'):
                interesting_subs.append(substrate)
    biolog_mappings = biolog_mappings.loc[interesting_subs, ]
    return biolog_mappings
    
    
    
def parallelizer(globber, function, nthreads):    
    
    globalpool = multiprocess.Pool(processes=nthreads, maxtasksperchild=1)
    
    all_files = list(glob.glob(globber))
    nitems_inchunk = int(len(all_files) / nthreads)
    if len(all_files) % nthreads !=0: nitems_inchunk += 1
    chunks = [all_files[x *nitems_inchunk : (x+1)* nitems_inchunk] for x in range(nthreads)]
    
    results = globalpool.imap(function, chunks, chunksize=1)
    
    results_df_all = []
    for result in results: 
        results_df_all.append(result)
    results_df_all = pnd.concat(results_df_all).reset_index(drop=True)
    
    globalpool.close()
    globalpool.join() 
    
    return results_df_all
        