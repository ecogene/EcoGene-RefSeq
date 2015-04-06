<?php 
class RefSeq{
	public function refseq_fna_file_read($file, &$sequence)
	{
		$fp = fopen($file, 'r');
		$head = fgets($fp);
		
		while (($line = fgets($fp))!== FALSE )
		{
			if($line[0]=='>')	continue;	
			$sequence .=  trim($line);
		}
		fclose($fp);
	}
	public function refseq_faa_file_read($file, $gi)
	{
		
		$fp = fopen($file, 'r');
		
		$header = "";
		$FASTA_seq = "";
		$flag = FALSE; // 
		while (($line = fgets($fp))!== FALSE )
		{
			if($line[0]=='>' && (strpos($line,$gi)!==FALSE))	{//start to read
				$header = $line;
				$flag = 1;		
				continue;		
			}elseif($line[0]=='>' && $flag) {  // stop to read when $flag = 2
				break;
			}
			
			if($flag==1){
				$FASTA_seq .= trim($line);
			}
		}
		fclose($fp);		
		 return array('main_string'=>$FASTA_seq, 'header'=>$header,);
		
		
	}
}
	
	