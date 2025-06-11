"""
# rotate
count!('total reads')
if read is [before:_ TTTACACTTTATGCTTCCGGCTCGTATGTTG after:_] => {
       	count!('rotated')
        rotated = before.concat(after)
       	# extract H and L
       	if rotated is [_ TTATTACTCGCTGCCCAACCAGCCATGGCC heavy:_ GCTTCCACCAAGGGCCCATCGGTCTTCCCG _] => {
               	count!('heavy')
               	if rotated is
                       	[_ CTGGCTGGTTTCGCTACGGCCGTGCAGGCA light:_ GGTCAGCCCAAGGCTGCCCCCTCGGTCACT _] => {
                               	count!('heavy + kappa')
                               	heavy.out!('heavy.fasta')
                               	light.out!('light.fasta')
                       	}
                       	[_ CTGGCTGGTTTCGCTACGGCCGTGCAGGCA light:_ GATCAAACGAAGGCTGCACCATCTGTCATT _] => {
                               	count!('heavy + lambda')
                               	heavy.out!('heavy.fasta')
                               	light.out!('light.fasta')
                       	}
       	}
}

# for the reverse complement
# rotate
if -read is [before:_ TTTACACTTTATGCTTCCGGCTCGTATGTTG after:_] => {
       	count!('rc rotated')
        rotated = before.concat(after)
       	# extract H and L
       	if rotated is [_ TTATTACTCGCTGCCCAACCAGCCATGGCC heavy:_ GCTTCCACCAAGGGCCCATCGGTCTTCCCG _] => {
               	count!('rc heavy')
               	if rotated is
                       	[_ CTGGCTGGTTTCGCTACGGCCGTGCAGGCA light:_ GGTCAGCCCAAGGCTGCCCCCTCGGTCACT _] => {
                               	count!('rc heavy + kappa')
                               	heavy.out!('heavy.fasta')
                               	light.out!('light.fasta')
                       	}
                       	[_ CTGGCTGGTTTCGCTACGGCCGTGCAGGCA light:_ GATCAAACGAAGGCTGCACCATCTGTCATT _] => {
                               	count!('rc heavy + lambda')
                               	heavy.out!('heavy.fasta')
                               	light.out!('light.fasta')
                       	}
       	}
}

"""