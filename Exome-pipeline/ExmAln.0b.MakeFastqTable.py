#!/home/local/users/jw/bin/python2.7
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Make FastqTable According to list of fastq names
#=========================================================================

from optparse import OptionParser
import re
# Take Fasetq File Name,
# Return RG: ID, SM, LB, PL, CN
name = re.compile('(OMG\d+-\d+)-') #PIPseq ID

#=========================================================================
# Modify this part each time run. Customized by different Batch of Data
#SM = '(OMG\d+-\d+-[A-Za-z0-9-]+)'
#ID = 1
#PL = 'Illumina'
#LB = 'Lib'
#LB_extra = '[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+' #160729_D00437_0328_BC9JMUANXX
#PAIR = 'Paired(1|2)'
#CN = 'OMG'
SM = '(CARE[a-zA-Z0-9]+)_'
ID = 1
PL = 'Illumina'
LB_extra = '(_DHE[a-zA-Z0-9]+)_'
#=========================================================================
SM = re.compile(SM)
#ID = re.compile(ID)
LB_extra = re.compile(LB_extra)
PAIR = re.compile(PAIR)

class FastqTableRecord:
    def __init__(self, f_path_1):
        self.Path1 = f_path_1
        self.Marker = self.Path1.rstrip('.gz').rstrip('.fq').rstrip('.fastq').rstrip('Paired1').rstrip('Paired2')
    def ParseName(self):
        self.Fname1 = self.Path1.split('\t')[-1]
        self.SM = SM.search(self.Fname1).group(1) # Sample Name
        self.ID = 1
        if 'Project_Clinical_WES' in self.Path1:
            self.LB = LB_extra.search(self.Path1).group(0)
        else:
            self.LB = LB
        self.CN = CN
        self.PL = PL
        self.RF1 = PAIR.search(self.Fname1).group(1) # Pair end #
    def checkRepeat(self, Sample_Dict):
        if self.SM in Sample_Dict:
            Sample_Dict[self.SM] += 1
            self.ID = Sample_Dict[self.SM]
        else:
            Sample_Dict[self.SM] = self.ID
    def search_mate(self, fq_list):
        for j in xrange(1, len(fq_list)):
            fq_path_2 = fq_list[j].strip()
            Marker = fq_path_2.rstrip('.gz').rstrip('.fq').rstrip('.fastq').rstrip('Paired1').rstrip('Paired2')
            if self.Marker == Marker:  # Find Mate
                print Marker
                self.Path2 = fq_path_2
                self.Fname2 = self.Path2.split('\t')[-1]
                self.RF2 = PAIR.search(self.Fname2).group(1)
                if self.RF1 == '1' and self.RF2 == '2':
                    self.FR, self.RE = self.Path1, self.Path2
                    fq_list.pop(j)
                    fq_list.pop(0)
                    return
                elif self.RF1 == '2' and self.RF2 == '1':
                    self.FR, self.RE = self.Path2, self.Path1
                    fq_list.pop(j)
                    fq_list.pop(0)
                    return 
                else:
                    print 'Error with File Name', fq_name, mate_name
                    exit()
        # Didn't find a mate:
        self.FR, self.RE = self.Path1, ""
        fq_list.pop(0)
        return
    def Format(self):
        RG = '{}\\t{}\\t{}\\t{}\\t{}\\t{}'.format('@RG', 'ID:' + str(self.ID), 'SM:' + self.SM, 'LB:' + self.LB, 'PL:' + self.PL, 'CN:' + self.CN)
        return '{}\t{}\t{}\n'.format(self.FR, RG, self.RE)

def MakeFastqTable(Input, Output):
    fq_list = [tmp.strip() for tmp in open(Input, 'rb').readlines()]
    fout = open(Output, 'wb')
    Sample_Dict = {}
    while len(fq_list) != 0:
        record = FastqTableRecord(fq_list[0].strip('\n'))
        record.ParseName()
        record.checkRepeat(Sample_Dict)
        record.search_mate(fq_list)
        fout.write(record.Format())
    fout.close()



def GetOptions():
    parser = OptionParser()
    parser.add_option('-i', '--input', dest='input',
                      metavar='input', help='Input fastq file listSRA')
    parser.add_option('-o', '--output', dest='output',
                      metavar='output', help='Output fastq table name')

    (options, args) = parser.parse_args()
    if options.output == None:
        options.output = 'Fastq_table.txt'
    return options.input, options.output

def main():
    Input, Output = GetOptions()
    MakeFastqTable(Input, Output)


if __name__ == '__main__':
    main()

def ParseName(fpath):
    fname = fpath.split('\t')[-1]
    _SM = SM.search(fname).group(1) # Sample Name
    _ID = ID.search(fname).group(1) # Read Group ID
    RF = PAIR.search(fname).group(1) # Pair end #
    return _ID, _SM, LB, PL, CN, RF
def MakePair_2(fq_list):
    while len(fq_list) > 0:
        fq_path = fq_list[0]
        fq_name = get_name(fq_path)
        R1, R2, j = search_mate(fq_name, fq_list)
        RG = get_RG(fq_list[0])
        if R2 != None:
            yield R1 + '\t' + RG + '\t' + R2 + '\n'
            fq_list.pop(j)
            fq_list.pop(0)
        else:
            yield R1 + '\t' + RG + '\t' + '\n'
def search_mate_2(fq_name, fq_list):
    fq_name = fq_name.rstrip('.fq.gz')
    fq_name = fq_name.split('_')

    for j in xrange(1, len(fq_list)):
        mate_name = get_name(fq_list[j]).rstrip('.fq.gz')
        mate_name = mate_name.split('_')
        flag = True
        for k in [0, 1, 2, 3]:
            if fq_name[k] != mate_name[k]:
                flag = False
                break
        if flag == True:  # Find Mate
            if fq_name[4] == '1' and mate_name[4] == '2':
                return fq_list[0], fq_list[j], j
            elif fq_name[4] == '2' and mate_name[4] == '1':
                return fq_list[j], fq_list[0], j
            else:
                print 'Error with File Name', '_'.join(fq_name), '_'.join(mate_name)
                exit()
    # Didn't find a mate:
    return '_'.join(fq_name), None, None
def MakePair(fq_list):
    count = 0
    while len(fq_list) > 0:
        fq_path = fq_list[0]
        fq_name = get_name(fq_path)
        R1, R2, j = search_mate(fq_name, fq_list)
        RG = get_RG(fq_list[0])
        if R2 != None:
            yield R1 + '\t' + RG + '\t' + R2 + '\n'
            fq_list.pop(j)
            fq_list.pop(0)
        else:
            yield R1 + '\t' + RG + '\t' + '\n'
            fq_list.pop(0)
def search_mate(fq_name, fq_list):
    fq_path_1 = "/".join(fq_list[0].split('/')[:-1])
    fq_name = fq_name.rstrip('.fastq.gz').rstrip('.fq.gz')
    _ID_1, _SM_1, LB_1, PL_1, CN_1 , RF_1= ParseName(fq_name)
    for j in xrange(1, len(fq_list)):
        fq_path_2 = "/".join(fq_list[j].split('/')[:-1])
        mate_name = get_name(fq_list[j]).rstrip('.fastq.gz').rstrip('.fq.gz')
        _ID_2, _SM_2, LB_2, PL_2, CN_2 , RF_2 = ParseName(mate_name)
        #print _ID_1, _SM_1, LB_1, PL_1, CN_1 , RF_1
        if fq_path_1 == fq_path_2 and _ID_1 == _ID_2 and _SM_1 == _SM_2 and LB_1 == LB_2 and PL_1 == PL_2 and CN_1 == CN_2:  # Find Mate
            print fq_path_1, _ID_1, _SM_1, LB_1, PL_1, CN_1 , RF_1
            print fq_path_2, _ID_2, _SM_2, LB_2, PL_2, CN_2 , RF_2
            if RF_1 == '1' and RF_2 == '2':
                return fq_list[0], fq_list[j], j
            elif RF_1 == '2' and RF_2 == '1':
                return fq_list[j], fq_list[0], j
            else:
                print 'Error with File Name', fq_name, mate_name
                exit()
    # Didn't find a mate:
    return fq_list[0], None, None
def get_RG(fq_name):
    fq_name = get_name(fq_name).split('.')[0]
    ID, SM, LB, PL, CN, RF = ParseName(fq_name)
    return '{}\\t{}\\t{}\\t{}\\t{}\\t{}'.format('@RG', 'ID:' + ID, 'SM:' + SM, 'LB:' + LB, 'PL:' + PL, 'CN:' + CN)
def MakeFastqTable_2(Input, Output):
    fq_list = [tmp.strip() for tmp in open(Input, 'rb').readlines()]
    fout = open(Output, 'wb')
    for record in MakePair(fq_list):
        fout.write(record)
    fout.close()
