import numpy as np
import wave
def txt_to_wav(in_file,out_file,sample_rate,sample_width,channels = 1):
    f = open(in_file,'rb')
    line = f.readline()
    data = []
    while line:   
        num = float(line.strip())  
        border = pow(2,sample_width * 8 - 1)
        data.append(max(min(border - 1,num*border),-(border-1)))  
        line = f.readline()
    f.close()
    wavdata = np.array(data)
    wavdata = wavdata.astype(np.short)

    
    f = wave.open(out_file, "wb")  
    f.setframerate(sample_rate)
    f.setsampwidth(sample_width)
    f.setnchannels(channels)
    f.writeframes(wavdata.tobytes())  # 将wav_data转换为二进制数据写入文件
    f.close()

path = 'xxx'
in_file = path + '\\xxx.txt'
out_file = path + '\\xxx.wav'
txt_to_wav(in_file,out_file,16000,2,1)
