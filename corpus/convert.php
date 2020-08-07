<?php
$source_dir=$argv[2];
$destination_dir=$argv[3];
$sample_rate=$argv[4];

echo "source=$source_dir\n";
echo "destination=$destination_dir\n";
echo "sample rate=$sample_rate\n";

$files=scandir ( $source_dir);
function endsWith($haystack, $needle)
{
    $length = strlen($needle);

    return $length === 0 || 
    (substr($haystack, -$length) === $needle);
}

for ($i=0;$i<sizeof($files);$i++){
	$source=$source_dir."/".$files[$i];
	$destination=$destination_dir."/".$files[$i];
	if (endsWith($source, ".wav")){
		$cmd="sox $source -c 1 -r $sample_rate $destination";
		echo "$cmd\n";
		$line=exec($cmd);
	}
}

?>
