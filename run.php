<?php

##################3
# args should be:
#    pattern = the pattern we should replace with #layers
#    min = minimal value of #layers
#    max = maximal value of #layers
#    command = the command to run
#    output = output filename to use
function binary_search($args) {
	 $low = $args['min'];
	 $high = $args['max'];
	  
	 while ($high >= $low) {
		$L = floor(($low + $high) / 2);

		$cmd = preg_replace("/{$args['pattern']}/", $L, $args['command']);

		$output = "{$args['output']}.$L";
		print("$cmd >& $output \n");
		system("$cmd > $output 2>&1");
		$ok = `grep 'test is OK'  $output `;

		if ($ok) {
			$high = $L - 1;
		} else {
			$low = $L + 1;
		}
	 }
	 return $L;
}



$min = 20;
foreach (Array(20, 40, 60, 80, 100) as $res) {
	$min = binary_search(
		Array(
			'min' => $min,
			'max' => 200,
			'command' => "./test_helib --res=$res --in=breast_cancer_classification.csv --t=12 --L=LLL ",
			'pattern' => 'LLL',
			"output" => "output.helib/output.helib.$res"
		)
	);
}

