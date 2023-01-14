import numpy as np
import wave
import glob
from pydub import AudioSegment as AS
import webrtcvad
import matplotlib.pyplot as plt
from scipy import signal

def vad_compute(byte_block,frame_size,sample_rate,sample_width,rank):
    vad = webrtcvad.Vad(rank)
    zeros_add = len(byte_block)%(frame_size*sample_width)
    for i in range(frame_size*sample_width - zeros_add):
        byte_block += b'\x00'
    
    frame_size = frame_size*sample_width # 帧bytes长
    labels = []
    mark = 1
    for pos in range(0, len(byte_block), frame_size):
        label = 1 if vad.is_speech(byte_block[pos:pos+frame_size], sample_rate=sample_rate) else 0
        if mark == 1:
            if label == 1:
                label = 0
            else:
                mark = 0
        labels.append(label)
    return (labels,byte_block)

def print_byte(out_file,byte_block,labels,is_pcm_wanted,is_vad_wanted):
    if is_pcm_wanted:
        with open(out_file + '\\aim.pcm', 'ab+') as pcmfile:
            pcmfile.write(byte_block)
            pcmfile.close()
    if is_vad_wanted:
        with open(out_file + '\\vad.txt', 'a+') as txtfile:
            a = ''
            for x in labels:
                a = a + str(x) +'\n' 
            txtfile.write(a)
            txtfile.close()

def audio_add(in_file,out_file,max_count,block_count,
                file_format = ['wav','flac','mp3'],
                sample_rate = 16000, sample_width = 2, frame_size = 160,
                get_byte = True, get_float = True):  
    
    is_pcm_wanted = True  
    is_vad_wanted = True
    is_draw_wanted = True
    is_print_wanted = False
    format_count = len(file_format)
    files = glob.glob(in_file + '\\**\\*.{0}'.format(file_format[0]),recursive=True)
    if format_count > 1:
        for i in range(format_count-1):
            files.extend(glob.glob(in_file + '\\**\\*.{0}'.format(file_format[i+1]), recursive=True))
            
    L = len(files)
    Lmax = min(L,max_count)
    byte_block = b''    
    float_block = []    
    for j in range(block_count):
        for i,f in enumerate(files[j*int(Lmax/block_count):(j+1)*int(Lmax/block_count)],start=1):
            print(f)
            if get_byte:
                bydata = AS.from_file(f,f.split(".")[1]).raw_data
                byte_block += bydata 
            if get_float:
                shdata = AS.from_file(f,f.split(".")[1]).get_array_of_samples()  
                fdata = np.array(shdata)/pow(2.0,sample_width * 8 - 1)
                float_block = np.concatenate([float_block,fdata])
        if is_vad_wanted & get_byte:
            (labels_1,zadd_byte_block) = vad_compute(byte_block,frame_size,sample_rate,sample_width,2)       
            (labels_2,_) = vad_compute(byte_block,frame_size,sample_rate,sample_width,3)
            labels = [0.5*(labels_1[k]+labels_2[k]) for k in range(len(labels_1))]
            
        if is_draw_wanted & get_float:
            plt.figure(dpi=100,figsize=(24,4))
            plt.plot(float_block)
            if is_vad_wanted & get_byte:
                Lab = np.array(labels)
                Lab = Lab.repeat(frame_size)
                plt.plot(Lab)                 
            else:
                labels = [];
            plt.show()
                      
        if is_print_wanted:
            if get_byte:
                if is_vad_wanted:
                    print_byte(out_file,zadd_byte_block,labels,is_pcm_wanted,is_vad_wanted)
                else:
                    print_byte(out_file,byte_block,labels,is_pcm_wanted,is_vad_wanted)
            if get_float:
                 print('\n')
                
        float_block = []
        byte_block = b''    
        print('{0}% accomplished'.format(round((j+1)/block_count*100,2)))
    return 0    
