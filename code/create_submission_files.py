import argparse
import os
import shutil

def copy_file_rec(file, out_dir, basedir=None, new_name=None):
    os.makedirs(out_dir, exist_ok=True)
    with open(file, 'r') as src:
        basefile = os.path.basename(file)
        if basedir is None:
            basedir = os.path.dirname(os.path.abspath(file))
        if new_name is None:
            new_name = os.path.join(out_dir, basefile)
        with open(new_name, 'w') as dest:
            for line in src.readlines():
                if line.find('\\input') == -1:
                    if line.find('\\includegraphics') == -1:
                        dest.write(line)
                    else:
                        fig_addr = line.split('{')[1].split('}')[0]
                        new_fig_addr = fig_addr.replace('/', '_')
                        shutil.copy(os.path.join(basedir, fig_addr), os.path.join(out_dir, new_fig_addr))
                        dest.write(line.replace(fig_addr, new_fig_addr))
                else:
                    addr = line.split('\\input{')[1].split('}')[0]
                    if not addr.endswith('.tex'):
                        addr += '.tex'
                    addr_path = os.path.join(basedir, addr)
                    new_addr = addr.replace('/', '_')
                    addr_name = os.path.join(out_dir, new_addr)
                    # print(addr_path, out_dir, addr_name, basedir)
                    copy_file_rec(addr_path, out_dir, new_name=addr_name, basedir=basedir)
                    dest.write(line.replace(addr.replace('.tex', ''), new_addr.replace('.tex', '')))


def create_submission_files(main_file, out_dir, extra_files=[]):
    os.makedirs(out_dir, exist_ok=True)
    for file in extra_files:
        shutil.copy(file, out_dir)
    copy_file_rec(main_file, out_dir)


def get_args():
    parser = argparse.ArgumentParser(description='Create flat latex directory for '
                                                 'submission')
    parser.add_argument('main_file', type=str,
                        help='Main .tex file.')
    parser.add_argument('out_dir', type=str,
                        help='Path to directory where the files should be generated.')
    parser.add_argument('-x', '--extra_files', nargs='+',
                        help='All additional files that should be copied.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    create_submission_files(args.main_file, args.out_dir, extra_files=args.extra_files)
