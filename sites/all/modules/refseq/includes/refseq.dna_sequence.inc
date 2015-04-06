<?php 
class DNASequence {
	
    protected $sequence;
    

 	public function __construct($sequence = NULL) {
    	if (empty($sequence)) {
      	throw new Exception('Cannot created an instance of DNASequence without a specified sequence.');
    	}
    	$this->sequence = strtoupper($sequence);
    	
  	}
  	
   /**
   * Returns the type of the entity.
   */
  	public function sequence() {    	
    	return $this->sequence;	    	
  	}
  
  /**
   * return the  reverse complement of the sequence
   * complement:  
   * A<->T 
   * C<->G 
   * B<->V 
   * D<->H 
   * K<->M 
   * Y<->R
   *
   */
  	
  	public function revc($string)
	{
			$string = strtoupper($string);

			$string = str_replace("Z", "T",(str_replace("T", "A",(str_replace("A", "Z", $string))))); /* A<->T */
			$string = str_replace("Z", "C",(str_replace("C", "G",(str_replace("G", "Z", $string))))); /* C<->G */
			
			$string = str_replace("Z", "V",(str_replace("V", "B",(str_replace("B", "Z", $string))))); /* B<->V */
			$string = str_replace("Z", "H",(str_replace("H", "D",(str_replace("D", "Z", $string))))); /* D<->H */
			
			$string = str_replace("Z", "M",(str_replace("M", "K",(str_replace("K", "Z", $string))))); /* K<->M */
			$string = str_replace("Z", "R",(str_replace("R", "Y",(str_replace("Y", "Z", $string))))); /* Y<->R */		
				
			$string = strrev($string);
			
			return $string;
	}    
	
	/* ambiguity code */

		/* B = C, G or T */

		/* D = A, G or T */

	 	/* H = A, C or T */

	 	/* K = G or T */

	 	/* M = A or C */

	 	/* N = A, C, G or T ------no change*/

	 	/* R = A or G (purines) */

	 	/* S = C or G ------no change**/

	 	/* V = A, C or G */

	 	/* W = A or T ------no change**/

	 	/* X = A, C, G or T */

	 	/* Y = C or T (pyrimidines) ------no change**/
	
	public function ambiguity_code(){
		
		$search = array('B', 'D', 'H', 'K', 'M', 'N' ,'R', 'S', 'V', 'W', 'X', 'Y');
		$replace   = array('CGT', 'AGT', 'ACT', 'GT', 'AC', 'ACGT' ,'AG', 'CG', 'ACG', 'AT', 'ACGT', 'CT');
		
		return array_combine($search, $replace);
		
	}
    
	
	public function gene_dna_sequence(&$sequence, $start, $end, $strand, $us=0,$ds=0)
	{
		$main_string = "";
		$up_string = "";
		$down_string = "";
		$seq_len = strlen($sequence);
		$us = intval($us);
		$ds = intval($ds);
		$main_string = substr($sequence, $start-1, $end- $start +1 );
		// clockwise gene
		
		if($strand=='-'){
			$us ^= $ds ^= $us ^= $ds;
			
			if(($start-$us)>0){
				
				$down_string = substr($sequence, $start - 1-$us, $us );
			
			}else{
				// circulal
				$down_string = substr($sequence, $start-1-$us). substr($sequence, 0, $start-1);
			}
			
			if(($end+$ds)>$seq_len){
					
				$up_string = substr($sequence, $end). substr($sequence, 0, $ds - ($seq_len - $end));
			}else{
					
				$up_string = substr($sequence, $end, $ds );
			}
			
			$main_string = DNASequence::revc($main_string);
			$up_string = DNASequence::revc($up_string);
			$down_string = DNASequence::revc($down_string);
			
		}else{
			
			if(($start-$us)>0){
				
				$up_string = substr($sequence, $start - 1 - $us, $us );
			
			}else{
				// circulal
				$up_string = substr($sequence, $start-1-$us). substr($sequence, 0, $start-1);
			}
			
			if(($end+$ds)>$seq_len){
					
				$down_string = substr($sequence, $end). substr($sequence, 0, $ds - ($seq_len - $end));
			}else{
					
				$down_string = substr($sequence, $end, $ds );
			}
			
		}
		
		return array('main_string'=>$main_string, 'up_string'=>$up_string, 'down_string'=>$down_string);
	
	}
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
	public function coor_sequence_format($sequence, $gene, $length=50, $triplet=10)
	{
		
		$FASTA_seq = "";
		$line_terminate = "\n";
		if(is_a($sequence,'string')){
			$header = "";
			$up_string = "";
			$down_string = "";
			$main_string  = $sequence;
			
		}else{			
			$up_string = isset($sequence['up_string'])? $sequence['up_string']:"";
			$down_string = isset($sequence['down_string'])? $sequence['down_string']:"";
			$main_string = isset($sequence['main_string'])? $sequence['main_string']:"";
		}
		// length of the up and down strings. 
		$uplen = strlen($up_string);
		$downlen = strlen($down_string);
		$mainlen = strlen($main_string);
		
		if ($gene['strand']=='+'){
			$up_string_coor = $gene['start']-strlen($up_string);
			$main_string_coor = $gene['start'];
			$down_string_coor = $gene['end']+1;
			$start = $gene['start'];
			$end = $gene['end'];
		}else{
			$up_string_coor = $gene['end'] + strlen($up_string);
			$main_string_coor = $gene['end'];
			$down_string_coor = $gene['start']-1;
			$start = $gene['end'];
			$end = $gene['start'];
		}
		$coor_orignal = "";
		$coor_relative = "";
		
		for ($i=1; $i<=$mainlen; $i+=$length){
			$t = (($i + $length)>$mainlen)? $mainlen : $i + $length-1;
			$coor_relative .= "[$i " ."$t]".$line_terminate;
			
			$t2 = ($gene['strand']=='+')? $start+$i-1 : $start-$i+1;
			$coor_orignal  .= "$t2".$line_terminate;
		}
		if($uplen>0){
			$coor_orignal = $line_terminate.$coor_orignal;
			$coor_relative = $line_terminate.$coor_relative;
			
			for ($i=1; $i<=$uplen; $i+=$length){
				if (($i + $length)>$uplen){
					$t = $uplen;
					$t2 = ($gene['strand']=='+')? $start-$uplen : $start+$uplen;
				}else{
					$t =  $i + $length-1;
					$t2 = ($gene['strand']=='+')? $start-$i-$length+1 : $start+$i+$length-1;
				}
				
				$coor_relative = "[-$t " ."-$i]".$line_terminate.$coor_relative;
				$coor_orignal  = "$t2".$line_terminate. $coor_orignal;
			}
		}
		
		if($downlen>0){
			$coor_orignal .= $line_terminate;
			$coor_relative .= $line_terminate;
			
			for ($i=1; $i<=$downlen; $i+=$length){
				$t = (($i + $length)>$downlen)? $downlen : $i + $length-1;
				$t1 = $i+$mainlen;
				$t += $mainlen;
				$coor_relative .= "[$t1 " ."$t]".$line_terminate;
			
				$t2 = ($gene['strand']=='+')? $end+$i : $end-$i;
				$coor_orignal  .= "$t2".$line_terminate;
			}
		}
		
		if(!empty($up_string))  {
			
			$up_string = str_pad($up_string, $uplen+(($length-$uplen%$length)%$length), ' ',STR_PAD_LEFT);
		}
		
		if(!empty($down_string))  {
			
			$down_string = str_pad($down_string, $downlen+(($length-$downlen%$length)%$length), ' ',STR_PAD_RIGHT);
			
		}			
		
		$str_seq = "";
		if(strlen($up_string)>0){
		for ($i=0; $i<strlen($up_string);$i++) {
	
			if (($i % $length)==0) { 
				$str_seq = $str_seq. $line_terminate;			
			
			}elseif (($triplet >0) && (($i % $triplet)==0)){

				$str_seq = $str_seq. " ";
			}
			
			$str_seq = $str_seq. $up_string[$i];
			
		}
				
		$str_seq .= $line_terminate;
		}
	
		
		for ($i=0; $i<strlen($main_string);$i++) {
	
			if (($i % $length)==0) { 
				
				$str_seq = $str_seq. $line_terminate;				
			
			}elseif (($triplet >0) && (($i % $triplet)==0)){

				$str_seq = $str_seq. " ";
			}			
			$str_seq = $str_seq. $main_string[$i];		
			
		}
				
		$str_seq .= $line_terminate;
		
		for ($i=0; $i<strlen($down_string);$i++) {
	
			if (($i % $length)==0) { 
								
				$str_seq = $str_seq. $line_terminate;			
				
			}elseif (($triplet >0) && (($i % $triplet)==0)){

				$str_seq = $str_seq. " ";
			}
			
			$str_seq = $str_seq. $down_string[$i];		
		}
				
		return array($coor_orignal, $str_seq, $coor_relative);
	}
	public function fasta_sequence_format($sequence, $length=75, $triplet=FALSE, $coor=FALSE)
	{
		
		$FASTA_seq = "";
		$line_terminate = "\n";
		if(is_a($sequence,'string')){
			$header = "";
			$up_string = "";
			$down_string = "";
			$main_string  = $sequence;
			
		}else{
			$header = isset($sequence['header'])? $sequence['header']:"";
			$up_string = isset($sequence['up_string'])? $sequence['up_string']:"";
			$down_string = isset($sequence['down_string'])? $sequence['down_string']:"";
			$main_string = isset($sequence['main_string'])? $sequence['main_string']:"";
		}
		
		if($triplet==FALSE){
				
			for ($i=0; $i<strlen($main_string);$i++) {
	
					if (($i % $length)==0) { 
	
						$FASTA_seq = $FASTA_seq. $line_terminate;
	
					}
					$FASTA_seq = $FASTA_seq. $main_string[$i];
			}
					
		}else{
				for ($i=0; $i<strlen($up_string);$i++) {
	
					if (($i % $length)==0) { 
	
					$FASTA_seq = $FASTA_seq.$line_terminate;
	
					}
					$FASTA_seq = $FASTA_seq. $up_string[$i];
				}
				$FASTA_seq = $FASTA_seq. $line_terminate;
				for ($i=0; $i<strlen($main_string);$i++) {
	
					if (($i % $length)==0) { 
	
					$FASTA_seq = $FASTA_seq. $line_terminate;
	
					}
					else
					{
						if (($triplet == true) && (($i % 3)==0)) {
	
							$FASTA_seq = $FASTA_seq. " ";
	
						}
					}
					$FASTA_seq = $FASTA_seq. $main_string[$i];
				}
				$FASTA_seq = $FASTA_seq. $line_terminate;
				for ($i=0; $i<strlen($down_string);$i++) {
	
					if (($i % $length)==0) { 
	
					$FASTA_seq = $FASTA_seq. $line_terminate;
	
					}
					$FASTA_seq = $FASTA_seq. $down_string[$i];
				}
			}
			
			if(empty($header))
				 return substr($FASTA_seq,1);
		    else 
		    	return $header.$line_terminate.substr($FASTA_seq,1);
	}

public function calculate_length($sequence,$type='aa'){

	if($type == "nt"){

		return strlen($sequence).' nt';

	} elseif ($type == "aa" ){

		return strlen($sequence).' aa';

	}elseif ($type == "other" ){

		return round(strlen($sequence)/3-1).' aa';

	}
	
}
public function calculateMW($sequence,$type='aa')
{
	

// Molecular Weight DNA
//Molecular Weight = (An x 329.21) + (Un x 306.17) + (Cn x 305.18) + (Gn x345.21) + 159.0
// Molecular Weight protein
// Protein Mw is calculated by the addition of average isotopic masses of amino acids in the protein and the average isotopic mass of one water molecule.
/*
monoisotopic average 
Alanine (A) 71.03711 71.0788 
Arginine (R) 156.10111 156.1875 
Asparagine (N) 114.04293 114.1038 
Aspartic acid (D) 115.02694 115.0886 
Cysteine (C) 103.00919 103.1388 
Glutamic acid (E) 129.04259 129.1155 
Glutae (Q) 128.05858 128.1307 
Glycine (G) 57.02146 57.0519 
Histidine (H) 137.05891 137.1411 
Isoleucine (I) 113.08406 113.1594 
Leucine (L) 113.08406 113.1594 
Lysine (K) 128.09496 128.1741 
Methionine (M) 131.04049 131.1926 
Phenylalanine (F) 147.06841 147.1766 
Proline (P) 97.05276 97.1167 
Serine (S) 87.03203 87.0782 
Threonine (T) 101.04768 101.1051 
Tryptophan (W) 186.07931 186.2132 
Tyrosine (Y) 163.06333 163.1760 
Valine (V) 99.06841 99.1326 
*/
	$arr_weight_rna = array(
				"A" => "329.21", 
				"U" => "306.17",
				"C" => "305.18",
				"G" => "345.21",	
				"O"  => "159.0"			
				);
	$arr_weight = array(

				"A" => "71.0788", 
				"R" => "156.1875",
				"N" => "114.1038",
				"D" => "115.0886",
				"C" => "103.1388",
				"E" => "129.1155",
				"Q" => "128.1307",
				"G" => "57.0519",
				"H" => "137.1411",
				"I" => "113.1594",
				"L" => "113.1594",
				"K" => "128.1741",
				"M" => "131.1926",
				"F" => "147.1766",
				"P" => "97.1167",
				"S" => "87.0782",
				"T" => "101.1051",
				"W" => "186.2132",
				"Y" => "163.1760",
				"V" => "99.1326",
				"H2O" => "18.01524"
				);			
	if($type=='aa')
	{		

		$protLen = strlen($sequence);
		$molecular_weight = $arr_weight["H2O"];

		for ($i=0; $i<$protLen;$i++)
		{
			$molecular_weight = $molecular_weight + $arr_weight[$sequence{$i}];
		}
	}
	
	if($type=='nt')
	{		
		
		$sequence = str_replace("T", "U", $sequence);
		
		$rnaLen = strlen($sequence);
		$molecular_weight = $arr_weight_rna["O"];

		for ($i=0; $i<$rnaLen;$i++)
		{
			$molecular_weight = $molecular_weight + $arr_weight_rna[$sequence{$i}];
		}	
		
	}

	return $molecular_weight;

}

}