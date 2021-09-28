from __future__ import division
import numpy as np
import operator
import math
import os
import io
import re
import random
import matplotlib.pyplot as plt
from random import shuffle

'''
A common experimental output in biomedical science is a list of genes implicated in a given biologicalprocess or disease. 
The results of a group of studies answering the same, or similar, questions can becombined by meta-analysis, aiming at 
finding a more reliable answer. 
Ranking aggregation methods canbe used for combining gene lists from various sources in meta-analyses. 
Evaluating a ranking aggregationmethod on a specific type of dataset before using it can support the reliability of the 
result since theproperty of a dataset can influence the performance of an algorithm. 
Evaluation of aggregation methodsis usually based on a simulated database because of the lack of a known truth for real data. 
Simulated datasets tend to be too small and neglect key features of experimental data, including heterogeneity ofquality, 
relevance and the inclusion of unranked lists. 
This script can generate simulated data to explore the performance of the aggregation methods as a function of emulating
the common scenarios of real genomics data, with various heterogeneity of quality, noise level, and a mix of unranked and 
ranked data using 20000 possible entities.  
It is based on the analysis of real genomic data and the simulated data generation model in the study of MAIC algorithm 
that samples a score from specific distribution for each entity to simulate the rankings. The investigated real data
include SARS-COV-2virus, cancer(NSCLC), and bacteria(macrophage apoptosis).

The model ranks entities by figure Z for each entity generated from Gaussian distribution.
For entity k in experiment(list) i,
Z_ki ~ N (mu_k, sigma_i^2)
log(sigma_i)~N(log(mean_noise_M), heterogeneity_D)
The length of each list is also sampled from specific distribution to emulating real datasets.

The generated data are lists and will be written in a file with name "(mean_noise)_(heterogeneity_D)_(dataset_type(
if unranked lists included)).txt"

Examples of running are shown at bottom of this file below "if __name__ == '__main__':"
'''


class SimulatedDataGenerator(object):
    """Class to generate simulated data"""

    def __init__(self, directory='./simulated_data/'):
        """
        list_number: Integer. the number of lists to generate in one input file

        mu_true: List(number of true entities) showing mean value of true entities
        mu_noise: List(number of true noise entities) showing mean value of noise entities

        directory: directory for generated file
        :param directory: directory for generated files

        """

        def exp_func(x, a, b, c):
            return a * np.exp(b * x) + c
        # m_true: number of true entities
        m_true = 1000
        # m_noise: number of true noise entities
        m_noise = 19000
        self.length_boundary_L = 2
        self.length_boundary_R = 20000
        # significance mu for each true entity, following a fitted curve.
        x = np.arange(1, m_true + 1)
        parameters = [0.1405463, -0.00990575, 0.06253942]
        y = [exp_func(x_k, *parameters) for x_k in x]
        ratio = 2.0 / y[0]
        self.mu_true = [(ratio * y_k) for y_k in y]
        self.mu_noise = np.zeros(m_noise)
        self.directory = directory
        # cutting_length: List(list_number) of integers to show the number of entities to left
        self.cutting_length = []
        self.is_ranked = []
        # [[number of ranked lists, number of unranked lists]]: [[11,0],[11,21],[5,0],[5,5]]
        # ranking_group[0] or large dataset and ranking_group[1] for small dataset
        self.ranking_group = [[11, 21], [5, 5]]
        # length parameter for real datasets
        self.sars_cov_2_length_mean_std = [[5.359720497751178, 1.805483212303656],
                                           [3.4364817135049583, 1.6937186948927982]]
        self.macrophage_NSCLC_length_mean_std = [[7.414072687219903, 1.8634137719951522],
                                                 [2.9905276829178913, 1.3389579490910093]]

        # heterogeneity of list quality and mean noise level.
        self.heterogeneity_D = [0.1, 0.5, 1, 3, 12]
        self.mean_noise_M = [0.5, 1, 3, 4, 12]
        self.D_plotting_mean_noise = [0.1, 0.5, 1, 3]
        self.M_plotting_D = [0.5, 1, 3, 4, 12]
        # create directory to store generated data
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

    def generate_lists(self, filenames, is_ranked, sigma, cutting_length, absent_rate=0):
        """
        generate simulated data in a file with certain parameters(for single running of the algorithm)
        :param filename: filename for generated file
        :param is_ranked: List(list_number) of boolean values to show if each generated list_i is ranked
        :param sigma: List(list_number) of numerical values to show experiment quality / noise level
        :return:
        """

        m_true = len(self.mu_true)
        m_noise = len(self.mu_noise)
        lists = []
        for i in range(len(is_ranked)):
            list_i = {}
            if sigma[i] == 0:
                for k in range(m_true):
                    list_i['entity' + str(k + 1)] = self.mu_true[k]
                noise_name = range(m_noise)
                for k in range(m_noise):
                    list_i['noise' + str(noise_name[k])] = self.mu_noise[k]
            else:
                for k in range(m_true):
                    list_i['entity' + str(k + 1)] = np.random.normal(loc=self.mu_true[k], scale=sigma[i])
                noise_name = range(m_noise)
                for k in range(m_noise):
                    list_i['noise' + str(noise_name[k])] = np.random.normal(loc=self.mu_noise[k], scale=sigma[i])
            ranked_items = sorted(list_i.items(), key=operator.itemgetter(1), reverse=True)
            ranked_list = [x[0] for x in ranked_items]
            # ranked_list = ranked_list[0:int(math.ceil(len(ranked_list) * self.cut_point_ratio[i]))]
            if (absent_rate > 0):
                # print('absent_rate:',absent_rate)
                gene_num_to_remove = self.length_boundary_R - int(cutting_length[i])
                gene_num_bottom = int(math.ceil(gene_num_to_remove * (1 - absent_rate)))
                gene_num_not_absent = int(cutting_length[i]) + gene_num_bottom
                gene_not_absent = random.sample(ranked_list, gene_num_not_absent)
                ranked_list = [x for x in ranked_list if x in gene_not_absent]
            ranked_list = ranked_list[0:int(cutting_length[i])]
            if is_ranked[i] != 1:
                shuffle(ranked_list)
            lists.append(ranked_list)
        # write generated data to ranked file
        f = open(filenames[0], "w")
        count = 1
        for i in range(len(lists)):
            if is_ranked[i] == 1:
                rank_mark = 'RANKED'
                f.write('\t'.join(('Category' + str(count), 'List' + str(count), rank_mark, 'NAMED_GENES',
                                   str(lists[i])[1:-1].replace(", ", "\t").replace("'", "") + '\n')))
            count += 1
        f.close()
        # write generated data to mix file
        f = open(filenames[1], "w")
        count = 1
        for i in range(len(lists)):
            if is_ranked[i] == 1:
                rank_mark = 'RANKED'
            else:
                rank_mark = 'UNRANKED'
            f.write('\t'.join(('Category' + str(count), 'List' + str(count), rank_mark, 'NAMED_GENES',
                               str(lists[i])[1:-1].replace(", ", "\t").replace("'", "") + '\n')))
            count += 1
        f.close()

    def log_gaussian(self, mean_noise, heterogeneity_d, ranking_group):
        """
        calculate sigmas by log(sigma_i)~N(log(mean_noise), h^2), heterogeneity_d = h^2
        :param mean_noise: mean noise (quality) for the generated data.
        :param heterogeneity_d: Heterogeneity of quality for generated data.
        :return:
        """
        ranked_num = ranking_group[0]
        unranked_num = ranking_group[1]
        list_num = ranked_num + unranked_num
        if mean_noise == 0:
            return np.zeros(list_num)
        log_means = math.log(mean_noise) * np.ones(list_num)
        # sd = heterogeneity_d * (abs(math.log(mean_noise))+1)
        sd = math.sqrt(heterogeneity_d)
        if sd != 0:
            noise_variance = np.random.normal(loc=0, scale=sd, size=list_num)
        else:
            noise_variance = np.zeros(list_num)
        sigmas = np.exp(log_means + noise_variance)
        return sigmas

    def log_norm_to_generate_list_length(self, mean_ln_length_mu, heterogeneity_ln_std, num):
        '''
        Generate length for each list independently given mean and std for the Normal distribution to generate ln value
        of C, which will convert to the list length after considering length boundaries..
        :param mean_ln_length_mu: mean value of Normal distribution in the generation process
        :param heterogeneity_ln_std: standard deviation used in the generation process
        :param num: number of lists
        :return: lengths
        '''
        list_num = num
        # log_means = math.log(mean_length_mu - self.length_boundary_L) * np.ones(list_num)
        log_means = mean_ln_length_mu * np.ones(list_num)
        sd = heterogeneity_ln_std
        if sd != 0:
            ln_variance = np.random.normal(loc=0, scale=sd, size=list_num)
        else:
            ln_variance = np.zeros(list_num)
        Cis = np.exp(log_means + ln_variance)
        lengths = np.floor(Cis + self.length_boundary_L)
        lengths[lengths > self.length_boundary_R] = self.length_boundary_R
        # print(Cis)
        # print(lengths)
        return lengths

    def length_generator(self, len_mean_std, ranking_group):
        '''
        Control the generation of length for simulated lists.
        :param len_mean_std: mean and standard deviation used in the process of generation
        :param ranking_group: ranking group related to dataset type and list numbers.
        :return: all_length
        '''
        ranked_num = ranking_group[0]
        unranked_num = ranking_group[1]
        ranked_length = self.log_norm_to_generate_list_length(len_mean_std[0][0], len_mean_std[0][1], ranked_num)
        if unranked_num > 0:
            unranked_length = self.log_norm_to_generate_list_length(len_mean_std[1][0], len_mean_std[1][1],
                                                                    unranked_num)
            all_length = np.append(ranked_length, unranked_length)
        else:
            all_length = ranked_length
        return all_length

    def preparation_and_generation(self, mean_noise, heterogeneity_d, ranking_group, additional_directory='',
                                   absent_rate=0):
        """
        Prepare for generation (the file name, sigmas, is_ranked) and generate data.
        :param additional_directory: the inner directory of file to generate.
        :param mean_noise: mean noise (quality) for the generated data.
        :param heterogeneity_d: Heterogeneity of quality for generated data.
        :param ranked_ratio: the ratio of data which is ranked, leaving the others as ranked data.
        :return: generated_file name
        """
        if absent_rate > 0:
            generated_file_rank = str(mean_noise) + '_' + str(heterogeneity_d) + '_' + str(ranking_group) + 'r_' + str(
                absent_rate) + '.txt'
            generated_file_mix = str(mean_noise) + '_' + str(heterogeneity_d) + '_' + str(ranking_group) + 'm_' + str(
                absent_rate) + '.txt'
        else:
            generated_file_rank = str(mean_noise) + '_' + str(heterogeneity_d) + '_' + str(ranking_group) + 'r.txt'
            generated_file_mix = str(mean_noise) + '_' + str(heterogeneity_d) + '_' + str(ranking_group) + 'm.txt'
        if len(self.is_ranked) > 0:
            is_ranked = self.is_ranked
        else:
            n_ranked = self.ranking_group[ranking_group][0]
            n_unrank = self.ranking_group[ranking_group][1]
            is_ranked = np.append(np.ones(n_ranked), np.zeros(n_unrank))
        if len(self.cutting_length) > 0:
            cutting_length = self.cutting_length
        else:
            if ranking_group == 0:
                len_mean_std = self.sars_cov_2_length_mean_std
            else:
                len_mean_std = self.macrophage_NSCLC_length_mean_std
            cutting_length = self.length_generator(len_mean_std, self.ranking_group[ranking_group])
        file_directory = self.directory + additional_directory
        if not os.path.exists(file_directory):
            os.makedirs(file_directory)
        filename_r = os.path.join(file_directory, generated_file_rank)
        filename_m = os.path.join(file_directory, generated_file_mix)
        sigmas = self.log_gaussian(mean_noise, heterogeneity_d, self.ranking_group[ranking_group])
        filenames = [filename_r, filename_m]
        self.generate_lists(filenames, is_ranked, sigmas, cutting_length, absent_rate=absent_rate)
        return filenames

    def generate_multiple_simulated_data(self, repeat_times=100, absent_rate=0):
        '''
        Generate multiple simulated data exploring various mean noise level and heterogeneity,
         sampling the each parameter setting for [repeat_times] times.
        :param repeat_times: repeated generation times for each setting of parameters.
        :param absent_rate: the rate for absent entities.
        :return:
        '''
        for i in range(repeat_times):
            print('Repeat:', i)
            for r in range(len(self.ranking_group)):
                # plotting the change of mean noise using different heterogeneity
                for d in self.D_plotting_mean_noise:
                    for m in self.mean_noise_M:
                        additional_file_directory = 'repeat_' + str(i + 1)
                        self.preparation_and_generation(m, d, r, additional_directory=additional_file_directory,
                                                        absent_rate=absent_rate)
                # plotting the change of heterogeneity using different mean noise
                for m in self.M_plotting_D:
                    for d in self.heterogeneity_D:
                        if d not in self.D_plotting_mean_noise:
                            additional_file_directory = 'repeat_' + str(i + 1)
                            self.preparation_and_generation(m, d, r, additional_directory=additional_file_directory,
                                                            absent_rate=absent_rate)

    def generate_multiple_simulated_data_with_absent_rates(self, repeat_times=100, absent_rates=[0.2, 0.5]):
        '''
        Generate multiple simulated data given absent rates
        :param repeat_times: repeated generation times for each setting of parameters.
        :param absent_rates: the rate for absent entities.
        :return:
        '''
        if len(repeat_times) == 1:
            for i in range(repeat_times):
                print('Repeat:', i)
                for r in range(len(self.ranking_group)):
                    for d in self.heterogeneity_D:
                        for absent_rate_i in absent_rates:
                            m = 3
                            additional_file_directory = 'repeat_' + str(i + 1)
                            self.preparation_and_generation(m, d, r, additional_directory=additional_file_directory,
                                                            absent_rate=absent_rate_i)
        elif len(repeat_times) == 2:
            for i in range(repeat_times[0], repeat_times[1]):
                print('Repeat:', i)
                for r in range(len(self.ranking_group)):
                    for d in self.heterogeneity_D:
                        for absent_rate_i in absent_rates:
                            m = 3
                            additional_file_directory = 'repeat_' + str(i + 1)
                            self.preparation_and_generation(m, d, r, additional_directory=additional_file_directory,
                                                            absent_rate=absent_rate_i)

    def cut_bottom_entities_by_highest_ranking_for_files(self, source_folder, result_folder, repeat_times=100):
        '''
        Cut bottom entities by highest ranking for files within the source_folder, searching files bu given parameters.
        :param source_folder: folder to store the data to cut
        :param result_folder: folder to store the result files
        :param repeat_times: repeated generation times of simulated data for each setting of parameters.
        :return:
        '''
        for i in range(repeat_times):
            print('Repeat:', i)
            for r in self.unranked_num:
                # plotting the change of mean noise using different heterogeneity
                for d in self.D_plotting_mean_noise:
                    for m in self.mean_noise_M:
                        additional_file_directory = 'repeat_' + str(i + 1)
                        file_name = str(m) + '_' + str(d) + '_' + str(r) + '.txt'
                        file_folder = os.path.join(source_folder, additional_file_directory)
                        new_file_path_folder = os.path.join(result_folder, additional_file_directory)
                        self.cut_bottom_entities_by_highest_ranking(file_folder, file_name, new_file_path_folder)
                # plotting the change of heterogeneity using different mean noise
                for m in self.M_plotting_D:
                    for d in self.heterogeneity_D:
                        if d not in self.D_plotting_mean_noise:
                            additional_file_directory = 'repeat_' + str(i + 1)
                            file_name = str(m) + '_' + str(d) + '_' + str(r) + '.txt'
                            file_folder = os.path.join(source_folder, additional_file_directory)
                            new_file_path_folder = os.path.join(result_folder, additional_file_directory)
                            self.cut_bottom_entities_by_highest_ranking(file_folder, file_name, new_file_path_folder)

    def cut_bottom_entities_by_highest_ranking(self, file_path_folder, file_name, new_file_path_folder, num_leanve=0
                                               , top_num=1000):
        """
        Remove entities not top ranked in any lists to leave only num_leanve entities if top_num =0.
        Remove entities not top ranked in any lists to leave only top entities with highest ranking<=top_num.
        :param file_path_folder: folder to store the file to process
        :param file_name: file name
        :param new_file_path_folder: folder to store result files
        :param num_leanve: leave only num_leanve entities if top_num =0
        :param top_num: leave only top entities with highest ranking<=top_num
        :return: length of ranked_entities left
        """
        file_path = os.path.join(file_path_folder, file_name)

        if not os.path.exists(new_file_path_folder):
            os.makedirs(new_file_path_folder)
        lists = self.read_file(file_path)
        if top_num > 0:
            ranked_entities = self.build_entity_ranks(lists, ranking_threshold=top_num)
        else:
            ranked_entities = self.build_entity_ranks(lists)
        if (num_leanve > 0):
            ranked_entities = ranked_entities[0:num_leanve]
        # write processed data to file
        new_file_path = os.path.join(new_file_path_folder, file_name)
        f = open(new_file_path, "w")
        for line in lists:
            columns = line.split("\t")
            # is_ranked = columns[2]
            # list_category = columns[0]
            # list_name = columns[1]
            list_to_write = columns[0:4]
            for entity_index in range(4, len(columns)):
                entity_name = columns[entity_index]
                if entity_name in ranked_entities:
                    list_to_write.append(entity_name)
            f.write('\t'.join(list_to_write) + '\n')
            # print(len(list_to_write), len(columns))
        f.close()
        return len(ranked_entities)

    def read_file(self, file_path):
        """
        Read the supplied file into memory
        Skips lines that begin with a hash character and stops reading the
        file when it encounters a line that begins with 5 or more hyphens
        """
        list_lines = []
        with io.open(file_path, "r") as input_file:
            for line in input_file.read().splitlines():
                if line.startswith("-----"):
                    break
                elif re.match("^\\s*$", line):
                    pass
                elif not line.startswith("#"):
                    list_lines.append(line)
        return list_lines

    def build_entity_ranks(self, list_lines, ranking_threshold=20000):
        '''
        Rank entities by the highest rank of each entity shown in a dataset.
        :param list_lines: lists
        :param ranking_threshold: only count entities ranked within the threshold
        :return:
        '''
        entity_dict = {}
        for line in list_lines:
            columns = line.split("\t")
            is_ranked = columns[2]
            for entity_index in range(4, len(columns)):
                entity_name = columns[entity_index]
                if is_ranked == 'RANKED':
                    entity_rank = entity_index - 3
                else:
                    entity_rank = len(columns) - 4
                if entity_name != '':
                    if entity_rank <= ranking_threshold:
                        if entity_name not in entity_dict:
                            entity_dict[entity_name] = entity_rank
                        elif entity_rank < entity_dict[entity_name]:
                            entity_dict[entity_name] = entity_rank
        ranked_list = []
        sorted_entities = sorted(entity_dict.items(), key=lambda item: item[1], reverse=False)
        print('original entity number:', len(sorted_entities))
        for entity in sorted_entities:
            ranked_list.append(entity[0])
            # if(len(ranked_list)==3000):
            #     print('3000th_highest_ranking:',entity)
        # print("entities number:",len(ranked_list))
        return ranked_list

    def cut_bottom_entities_by_highest_ranking_for_files_in_folders(self, source_folder, result_folder,
                                                                    repeat_times=100):
        '''
        Cut bottom entities by highest ranking for files in folders, by searching files within given folder.
        :param source_folder: source folder including /repeat_[i]/[].txt files.
        :param result_folder: result folder
        :param repeat_times: repeated generation times of simulated data for each setting of parameters.
        :return:
        '''
        file_names = []
        file_name_read = False
        number_of_files = 0
        number_of_entities_left = []
        for i in range(repeat_times):
            print('Repeat:', i)
            additional_file_directory = 'repeat_' + str(i + 1)
            file_folder = os.path.join(source_folder, additional_file_directory)
            new_file_path_folder = os.path.join(result_folder, additional_file_directory)
            if not file_name_read:
                file_names = [x for x in os.listdir(file_folder) if ((not x.startswith('.')) and ('m' not in x))]
                file_names.sort()
                file_name_read = True
            for file_name in file_names:
                entities_num_left = self.cut_bottom_entities_by_highest_ranking(file_folder, file_name,
                                                                                new_file_path_folder)
                number_of_files = number_of_files + 1
                number_of_entities_left.append(entities_num_left)
        print('number_of_files:', number_of_files)
        print(number_of_entities_left)
        print('largest entities number:', max(number_of_entities_left))
        # plot number of entities left
        fig, ax = plt.subplots(1, 1)
        ax.hist(number_of_entities_left, bins=100, density=False)
        ax.set_xlabel('number of genes')
        ax.set_ylabel('number of files')
        plt.title('number of genes left after cutting not top 1000 ranked genes')
        ax.legend()
        plt.savefig(os.path.join(result_folder, 'gene_number_left_plotting.pdf'), format='pdf')
        plt.show()

    def maic_2_bard_topk_format(self, file_path_folder, file_name, new_file_path_folder):
        '''
        Transfer MAIC formated input file to BARD required format
        :param file_path_folder: the folder including the file to process
        :param file_name: file name
        :param new_file_path_folder: result folder
        :return:
        '''
        file_path = os.path.join(file_path_folder, file_name)
        if not os.path.exists(new_file_path_folder):
            os.makedirs(new_file_path_folder)
        lists = self.read_file(file_path)
        # write processed data to file
        new_file_path = os.path.join(new_file_path_folder, file_name)
        f = open(new_file_path, "w")
        all_entities = []
        ranker_names = []
        all_entities_len = []
        max_len = 0
        for line in lists:
            entities = line.split("\t")
            ranker_names.append(entities[1])
            genes = entities[4:]
            all_entities.append(genes)
            len_x = len(genes)
            all_entities_len.append(len_x)
            if max_len < len_x:
                max_len = len_x
        list_number = len(all_entities)
        f.write('Rank\t' + '\t'.join(ranker_names) + '\n')
        for index in range(0, max_len):
            string_to_write = str(index + 1) + '\t'
            for j in range(list_number):
                if index < all_entities_len[j]:
                    string_to_write = string_to_write + all_entities[j][index]
                if j < list_number - 1:
                    string_to_write = string_to_write + '\t'
            f.write(string_to_write + '\n')
        f.close()
        return

    def maic_2_bard_topk_format_in_folders(self, source_folder, result_folder, repeat_times=100):
        '''
        transfer MAIC formated input file to BARD required format for multiple files within source_folder
        :param source_folder: the source filder including files to process
        :param result_folder: result folder to store processed data
        :param repeat_times: repeated generation times of simulated data for each setting of parameters.
        :return:
        '''
        file_names = []
        file_name_read = False
        number_of_files = 0
        for i in range(repeat_times):
            print('Repeat:', i)
            additional_file_directory = 'repeat_' + str(i + 1)
            file_folder = os.path.join(source_folder, additional_file_directory)
            new_file_path_folder = os.path.join(result_folder, additional_file_directory)
            if not file_name_read:
                file_names = [x for x in os.listdir(file_folder) if ((not x.startswith('.')) and ('m' not in x))]
                file_names.sort()
                file_name_read = True
            for file_name in file_names:
                self.maic_2_bard_topk_format(file_folder, file_name, new_file_path_folder)
                number_of_files = number_of_files + 1
        print('number_of_files:', number_of_files)


# example to run
if __name__ == '__main__':
    generator = SimulatedDataGenerator(directory='./simulated_data/')
    # -generate simulated data with various mean noise level and heterogeneity.
    generator.generate_multiple_simulated_data(repeat_times=3)
    # -generate simulated data with various absent rate
    generator.generate_multiple_simulated_data_with_absent_rates(repeat_times=[94, 100])
