import os
import platform
SYS = platform.system()

def merge_reads_command(reads_dir, filetype, target_dir_path, file_name):
    # This already assumes the read directory and target directory are valid
    cmd = ""
    subcommand = ""
    os_command = ""
    output_filename = "{}{}.pod5".format(target_dir_path, file_name)
    
    if filetype == 'fast5':
        subcommand = 'convert fast5'
    elif filetype == 'pod5':
        subcommand = ''
    
    if SYS == 'Windows':
        os_command = '--recursive {}'.format(reads_dir)
    else:
        os_command = '{}*.{}'.format(reads_dir, filetype)
        
    cmd = "pod5 {} --force-overwrite {} -o {}".format(subcommand, os_command, output_filename)

    return cmd


def validate_read_directory(reads_dir):
    directory_exists = os.path.isdir(reads_dir)
    homogenous_files = True
    filetype = ""
    if directory_exists:
        directory_files = os.listdir(reads_dir)
        ext_list = []
        for file in directory_files:
            ext = file.split('.')[-1]
            ext_list.append(ext)
        uniques = list(set(ext_list))
        if (len(uniques) != 1):
            homogenous_files = False
            print('Passed reads directory not homogenous. Filetypes found: {}'.format(uniques))
        else:
            filetype = uniques[0]
            
    return filetype


def validate_target_directory(target_dir):
    directory_exists = os.path.isdir(target_dir)
    return directory_exists
        
    
def generate_merged_pod5(reads_dir, target_dir_path, file_name):
    filetype = validate_read_directory(reads_dir)
    valid_target = validate_target_directory(target_dir_path)
    st = 1
    if (valid_target and (filetype != "")):
        cmd = merge_reads_command(reads_dir, filetype, target_dir_path, file_name)
        st = os.system(cmd)
    return (st, filetype, valid_target)
    