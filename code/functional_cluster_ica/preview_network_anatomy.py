#!/usr/bin/env python3
"""Review published fixed functional ICA anatomy/task contrasts without 4D data."""
import argparse
import os
from pathlib import Path
import network_summary as s


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sm0-record',type=Path,required=True);p.add_argument('--sm5-record',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--fsldir',type=Path,default=Path(os.environ.get('FSLDIR','/usr/local/fsl')));args=p.parse_args()
    sources,states,rows,ref,mask,original,families=s.load_sources(args.sm0_record,args.sm5_record)
    args.output.mkdir(parents=True,exist_ok=True);(args.output/'figures').mkdir(exist_ok=True)
    s.anatomy(families,ref,mask,args.output,args.fsldir);s.task_summary(sources,args.output)
    s.f.save(args.output/'preview_provenance.json',dict(source_records={k:str(v.resolve()) for k,v in sources.items()},component_maps={k:v['map_sha256'] for k,v in families.items()},anatomy_components=32,participant_4D_data_read=False,new_regression_fits=False,task_results='Algebraic paired contrasts of previously standardized participant estimates'))
    print('Anatomy and paired task contrasts ready:',args.output)


if __name__=='__main__':main()
