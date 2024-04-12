import os
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import xf_params as xfp

def dorado_command(basecaller_path, model_type, reference_path, pod5_path, out_path, out_name):
    """
    dorado_command generates a command to run the basecaller dorado.
    This instance of dorado will perform both basecalling and 
    alignment. 
    
    Parameters:
    basecaller_path: path to the basecaller, as str
    model_type: either one of the auto selectd models when 'fast', 'hac', or 'sup' inputted or a specific model, as str
    reference_path: file pathway to consensus or predicted ground truth reference fasta, as str
    pod5_path: path to the pod5 file to be basecalled, as str
    out_path: path to the output directory, as str
    out_name: name the basecalled fq file will be given, as str
    
    Returns:
    a command string.
    """
    cmd = "{} basecaller {} --no-trim --emit-moves --secondary n --reference {} {} > {}{}.bam".format(basecaller_path, model_type, reference_path, pod5_path, out_path, out_name)

    print('XenoFind [STATUS] - basecalling command generated: "{}"'.format(cmd))
    return cmd
