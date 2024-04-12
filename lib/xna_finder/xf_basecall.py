import os
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import xf_params as xfp

def dorado_command(basecaller_path, model_type, reference_path, pod5_path, out_path, out_name):
    """
    dorado_command generates a command to run the basecaller dorado.
    
    Parameters:
    basecaller_path: path to the basecaller, as str
    pod5_path: path to the pod5 file to be basecalled, as str
    out_path: path to the output directory, as str
    out_name: name the basecalled fq file will be given, as str
    
    Returns:
    a command string.
    """
    cmd = "{} basecaller {} --no-trim --emit-moves --reference {} {} > {}{}.fq".format(basecaller_path, model_type, reference_path, pod5_path, out_path, out_name)
    bam_path = out_path+out_name
    '''
    need to add removal of secondary and unalign read flags to this command 
    '''
    print('[Basecalling]: Command Generated: "{}"'.format(cmd))
    return cmd, bam_path
