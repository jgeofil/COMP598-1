import sys, re, math
# ============================================ Student info methods================================================
def get_student_name():
	# @TO_STUDENT: Write your name here
	student_name = "Jeremy Georges-Filteau"
	if not student_name:
		raise ValueError("Error: you forgot to add your name in get_student_name method.")
	return student_name

def get_student_id():
	# @TO_STUDENT: Write your student id here
	student_id = "260713547"
	if not student_id:
		raise ValueError("Error: you forgot to add your student id in get_student_id method.")
	return student_id
# =================================================================================================================

# =================================== Validate input/output methods================================================
def validate_Q3_1_input_format(subopt_result):
	if not isinstance(subopt_result, list) or [sr for sr in subopt_result if re.search("[^\(\)\.]", sr)]:
		raise ValueError("Error: your input should be list of strings (secondary structures, ie alphabet is '().')")

def validate_Q3_1_output_format(result):
	if not isinstance(result, list) or [sr for sr in result if not isinstance(sr, list)]:
		raise ValueError("Error: your output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ].")
	if [sr for sr in result if sr[0] >= sr[1]]:
		raise ValueError("Error: check your i and j values. They should satisfy: i > j.")
# =================================================================================================================
# ================================================== Helper methods================================================
def parse_subopt_result_file(filepath):
	'''
	Parsing of a standard txt file that contains result of
		"RNAsubopt -p __k__ < myFasta.fasta > subopt_result_file.txt"
		where __k__ is parameter. (Filename chosen randomly. Please, feel free to use your own names.)
	@args:
	filepath: (full or relative) path to the subopt_result_file.txt.
	(Filename chosen randomly. Please, feel free to use your own names.)

	@return: subopt_result: list of the strings (assumed to be secondary structures)
	'''
	subopt_result = []
	with open(filepath, 'r') as f:
		for i, line in enumerate(f):
			if i < 2:
				continue
			subopt_result.append(line.strip())
	return subopt_result

def parse_dot_ps_file(filepath):
	'''
	Parsing of a dot.ps file that contains result of RNAfold program
	@args:
	filepath: (full or relative) path to the dot.ps.
	@return:
	dot_ps_result: list f lists with i, j, freq_i_j
	'''
	dot_ps_result = []
	with open(filepath, 'r') as f:
		is_data = False
		for line in f:
			if not is_data and line.startswith('%start of base pair probability data'):
				is_data = True
				continue
			elif is_data and line.startswith('showpage'):
				break
			elif is_data:
				# take only first 3 numbers
				data_line = line.split()[:3]
				dot_ps_result.append(
					[int(data_line[0]), int(data_line[1]), float(data_line[2])]
				)
	return dot_ps_result

# =================================================================================================================
def get_answer_Q3_1(subopt_result):
	'''
	This method should be implemented by student.
	@args:
	subopt_result: a list of the secondary structures produced by RNAsubopt -p <k> for particular input

	@return:
	result: list of lists (example is below) with indices and relevant frequency.
	example [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16], ....... ]

	@note:
	check input/output as advised in code. Question will be marked as 0 in case of not following the formats.
	'''
	# basic check for the proper input
	validate_Q3_1_input_format(subopt_result)
	# @TO_STUDENT: Write your code here

	seqs = len(subopt_result)
	seqLen = len(subopt_result[0])
	openings = [[] for x in range(0, seqs)]
	occurences = [[0]*seqLen for x in range(0, seqLen)]

	for n in range(0, seqLen):
		for m in range(0, seqs):
			symbol = subopt_result[m][n]
			if symbol == '(':
				openings[m].append(n)
			elif symbol == ')':
				occurences[openings[m].pop()][n] += 1
	result = []
	for m, narr in enumerate(occurences):
		for n, val in enumerate(narr):
			if val != 0:
				result.append([m,n,val/float(seqs)])

	# @TO_STUDENT: use result variable for results. below is an example of an expected format for result.
	#result = [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16] ]

	# @TO_STUDENT: output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ]
	# use validate_Q3_output_format(result) to validate the output
	validate_Q3_1_output_format(result)
	return result

def get_answer_Q3_2(q3_1_result, dot_ps_result):
	'''
	This method should be implemented by student.
	Compare output from RNAfold and result of question3_1 for the same sequence and return an error (see text assignment)
	result_error is expected to be numeric
	'''
	result_error = 0
	# @TO_STUDENT: Write your code here (trust me, answer is not 0 :-) )

	def calcError(a,b):
		return (a-b)**2

	idx1, idx2 = 0,0
	end1, end2 = len(q3_1_result), len(dot_ps_result)

	while idx1 < end1 or idx2 < end2:
		if idx1 >= end1:
			val2 = dot_ps_result[idx2]
			result_error += calcError(0, val2[2])
			idx2 += 1
		elif idx2 >= end2:
			val1 = q3_1_result[idx1]
			result_error += calcError(val1[2], 0)
			idx1 += 1
		else:
			val1, val2 = q3_1_result[idx1], dot_ps_result[idx2]
			if val1[0] == val2[0] and val1[1] == val2[1]:
				result_error += calcError(val1[2], val2[2])
				idx1 += 1
				idx2 += 1
			elif val1[0] == val2[0]:
				if val1[1] > val2[1]:
					result_error += calcError(0, val2[2])
					idx2 += 1
				else:
					result_error += calcError(val1[2], 0)
					idx1 += 1
			elif val1[0] > val2[0]:
				result_error += calcError(0, val2[2])
				idx2 += 1
			elif val1[0] < val2[0]:
				result_error += calcError(val1[2], 0)
				idx1 += 1
	result_error = math.sqrt(result_error)
	print result_error
	return result_error

# @TO_STUDENT: You can test your methods below by calling methods. Workflow is given already (can be changed).
# @TO_STUDENT: Everything below this point will not be considered as a solution and will be deleted for grading.
# @TO_STUDENT: Advise: be careful with idents. Use only tabs, or only FOUR spaces. NEVER mix them.

print("This is a solution of %s, student_id is %s" % (get_student_name(), get_student_id()) )

subopt_result_filepath = "out.txt"
dot_ps_filepath = "dot.ps"

# parsing RNAsubopt result file
subopt_result = parse_subopt_result_file(subopt_result_filepath)

# solving quesion Q3_1
q3_1_result = get_answer_Q3_1(subopt_result)

# parsing dot.ps file
dot_ps_result = parse_dot_ps_file(dot_ps_filepath)

# solving question Q3_2
q3_2_result = get_answer_Q3_2(q3_1_result, dot_ps_result)
