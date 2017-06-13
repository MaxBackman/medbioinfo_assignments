
my_file = open("input/blah.fq", "r")
file_contents = my_file.read()
# print(file_contents)

line_array = file_contents.split('\n')

print(int(len(line_array)/4))


for i in list(range(1, int(len(line_array)/4)+1)): # creates an array of 1,2,3... to the number of reads in the FastQ.fa file
    utan_snabela = line_array[i*4-4][1:len(line_array[i*4-4])] # remove initial '@'
    print(utan_snabela)

my_file.close()
