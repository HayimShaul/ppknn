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
		$ok = `grep 'test is ok'  $output `;

		if ($ok) {
			$high = $L - 1;
		} else {
			$low = $L + 1;
		}
	 }
	 return $L;
}

#binary_search(
#	Array(
#		'min' => 85,
#		'max' => 110,
#		'command' => "./test_helib --t=12 --L=LLL --in=boston_hotels.json",
#		'pattern' => 'LLL',
#		"output" => "output.helib"
#	)
#);

$settings = Array(
#		Array("n" => 500, "p" => 70, "L" => 71),
#		Array("n" => 500, "p" => 80, "L" => 93),
#		Array("n" => 500, "p" => 90, "L" => 80),
#		Array("n" => 500, "p" => 100, "L" => 80)
	);

#foreach ($settings as $s) {
#	$n = $s['n'];
#	$p = $s['p'];
#	binary_search(
#		Array(
#			'min' => $s['L'],
#			'max' => $s['L'],
#			'command' => "./test_helib --n=$n --p=$p --t=12 --L=LLL ",
#			'pattern' => 'LLL',
#			"output" => "output.helib.$n.$p"
#		)
#	);
#}


for ($n = 1000; $n < 10000; $n += 1000) {
#	for ($p = 50; $p < 110; $p += 10) {
	foreach (Array(70,90,100)  as $p) {
		$found = false;

		foreach ($settings as $s) {
			if (($n == $s['n']) && ($p == $s['p']))
				$found = true;
		}

		if (!$found) {
			binary_search(
				Array(
					'min' => 70,
					'max' => 120,
					'command' => "./test_helib --n=$n --p=$p --t=12 --L=LLL ",
					'pattern' => 'LLL',
					"output" => "output.helib.$n.$p"
				)
			);
		}
	}
}

#for ($n = 500; $n < 1000; $n += 100) {
#	for ($p = 50; $p < 110; $p += 10) {
#		$found = false;
#
#		foreach ($settings as $s) {
#			if (($n == $s['n']) && ($p == $s['p']))
#				$found = true;
#		}
#
#		if (!$found) {
#			binary_search(
#				Array(
#					'min' => 70,
#					'max' => 110,
#					'command' => "./test_helib --n=$n --p=$p --t=12 --L=LLL ",
#					'pattern' => 'LLL',
#					"output" => "output.helib.$n.$p"
#				)
#			);
#		}
#	}
#}
