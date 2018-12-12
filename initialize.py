import os

bowtie2_path = '/catt/bowtie2-2.3.4.3-linux-x86_64/bowtie2'
bowtie2_build_path = '/catt/bowtie2-2.3.4.3-linux-x86_64/bowtie2-build'
samtool_path = '/catt/samtools-1.9/samtools'

if __name__ == '__main__':

    refs = os.listdir("./resource")
    for name in refs:
        if name[0]=='T':
            os.system("%s ./resource/%s ./resource/%s" % (bowtie2_build_path, name, name))

    os.system("gcc ./resource/hamming.c -shared -o ./resource/hamming.so")
