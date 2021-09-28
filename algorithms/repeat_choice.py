import io
import random
import re
import json
import numpy as np

class Repeat_Choice(object):
    """Class to manage the repeat choice method for ranking aggregation
    """

    def __init__(self, file_path, output_file='none'):
        self.entity_number = 0
        self.list_number = 0
        self.list_lines = []
        self.entities = []
        self.read_file(file_path)
        self.initialize_parameters()
        self.output_file = output_file
        self.ranked_entities = []

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

    def initialize_parameters(self):
        list_index = 0
        self.list_number = len(self.list_lines)
        for line in self.list_lines:
            columns = line.split("\t")
            entity_number_in_list = len(columns) - 4
            for entity_index in range(4, len(columns)):
                rank = entity_index - 3
                entity_name = columns[entity_index]
                if entity_name != '':
                    if entity_name not in self.entities:
                        self.entities.append(entity_name)
            list_index = list_index + 1
        self.entity_number = len(self.entities)

    def run_repeat_choice(self):
        '''
        Run the algorithm
        :return:
        '''
        self.ranked_entities = []
        self.ranked_entities.append(self.entities)
        list_index_to_choice = list(range(self.list_number))
        while len(self.ranked_entities) < self.entity_number and len(list_index_to_choice)>0:
            k = np.random.randint(len(list_index_to_choice))
            index_k = list_index_to_choice[k]
            chosen_list = self.list_lines[index_k]
            list_index_to_choice.remove(index_k)
            self.ranked_entities = self.merge_order(self.ranked_entities, chosen_list)
        self.ranked_entities = self.merge_order(self.ranked_entities, separated=True)

    def merge_order(self, ordered_entities, line = [], separated = False):
        merged_order = []
        if separated:
            for group in ordered_entities:
                random.shuffle(group)
                for e in group:
                    merged_order.append(e)
            return merged_order
        columns = line.split("\t")
        merged_order = []
        is_ranked = columns[2]
        entities_in_line = columns[4:]
        for group in ordered_entities:
            entities_left = group
            entities_ordered = []
            if is_ranked == 'RANKED':
                for e in entities_in_line:
                    if e in group:
                        merged_order.append([e])
                        entities_left.remove(e)
                merged_order.append(entities_left)
            else:
                for e in group:
                    if e in entities_in_line:
                        entities_ordered.append(e)
                        entities_left.remove(e)
                merged_order.append(entities_ordered)
                merged_order.append(entities_left)
        return merged_order

    def output_result(self):
        if self.output_file != 'none':
            try:
                result_file = open(self.output_file + '.txt', 'w')
                result_file.write('\t'.join(self.ranked_entities))
                result_file.close()
            except Exception as e:
                print('error when writting the result of repeat choice: '+str(e))
        return self.ranked_entities


if __name__ == '__main__':
    repeat_choice = Repeat_Choice('/Users/s1718825/Documents/MAIC/ExistingMethods/1_3_1.txt',output_file='/Users/s1718825/Documents/MAIC/ExistingMethods/RC_result.txt')
    repeat_choice.run_repeat_choice()
    result = repeat_choice.output_result()
    print(result)
