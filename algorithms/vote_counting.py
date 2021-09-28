import io
import re
import json


class Vote_counting(object):
    """Class to manage the vote counting method"""

    def __init__(self, file_path, output_file='none'):
        self.list_lines = []
        self.entity_dict = {}
        self.read_file(file_path)
        self.build_entity_dic()
        self.output_file = output_file

    def read_file(self, file_path):
        """
        Read the supplied file into memory
        Skips lines that begin with a hash character and stops reading the
        file when it encounters a line that begins with 5 or more hyphens
        """
        with io.open(file_path, "r") as input_file:
            for line in input_file.read().splitlines():
                if line.startswith("-----"):
                    break
                elif re.match("^\\s*$", line):
                    pass
                elif not line.startswith("#"):
                    self.list_lines.append(line)

    def build_entity_dic(self):
        '''
        Build entity dictionary for appearing times.
        :return:
        '''
        list_of_lists = []
        for line in self.list_lines:
            columns = line.split("\t")
            for entity_index in range(4, len(columns)):
                entity_name = columns[entity_index]
                if entity_name != '':
                    if entity_name not in self.entity_dict:
                        self.entity_dict[entity_name] = 1
                    else:
                        self.entity_dict[entity_name] = self.entity_dict[entity_name] + 1

    def run_vc(self):
        '''
        Run the algorithm
        :return:
        '''
        ranked_list = []
        sorted_entities = sorted(self.entity_dict.items(), key=lambda item: item[1], reverse=True)
        for entity in sorted_entities:
            ranked_list.append(entity[0])
        if self.output_file != 'none':
            try:
                fileObject = open((self.output_file + 'frequency.txt'), 'w')
                fileObject.write(json.dumps(sorted_entities))
                fileObject.close()
                result_file = open(self.output_file + '.txt', 'w')
                result_file.write('\t'.join(ranked_list))
                result_file.close()
            except Exception as e:
                print('error when writting the result of vote counting: '+str(e))
        return ranked_list


if __name__ == '__main__':
    vc = Vote_counting(
        'D:\programmer\PythonWorkSpace\crossvalidation_develop\\real_data\host.noCRISPR20-reformated.txt',
        output_file='D:\programmer\PythonWorkSpace\crossvalidation_develop\\real_data\\vc_result')
    result = vc.run_vc()
