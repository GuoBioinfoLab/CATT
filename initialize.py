import os

bowtie2_path = '/Users/kroaity/定期删除/bowtie2-2.3.4.3-macos-x86_64/bowtie2'
bowtie2_build_path = '/Users/kroaity/定期删除/bowtie2-2.3.4.3-macos-x86_64/bowtie2-build'
samtool_path = 'samtools'

if __name__ == '__main__':

    refs = os.listdir("./resource")
    for name in refs:
        os.system("%s ./resource/%s ./resource/%s" % (bowtie2_build_path, name, name))

    os.system("gcc ./resource/hamming.c -shared -o ./resource/hamming.so")
