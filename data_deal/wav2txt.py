import numpy as np
from pydub import AudioSegment as AS

def wav_to_txt(in_file,out_file,sample_width):
    shdata = AS.from_file(in_file,'{0}'.format(in_file.split('.')[1])).get_array_of_samples()
    str_data = ''
    for i in range(len(shdata)):
        str_data += str('{:f}'.format(float(shdata[i])/pow(2,sample_width * 8 - 1))) + '\n'
    with open(out_file, 'w') as txtfile:
        txtfile.write(str_data)
        txtfile.close()
            
path = 'xxx'
in_file = path + '\\xxx.wav'
out_file = path + '\\xxx.txt'
wav_to_txt(in_file,out_file,2)
