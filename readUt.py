def clean_and_split(line):
	# ignore anything that follows after a '#' as a comment 
	line        = line.partition('#')[0]
	
	# split the line in Arguments
	line        = line.rstrip()
	lineArgs    = line.split()

	return lineArgs
