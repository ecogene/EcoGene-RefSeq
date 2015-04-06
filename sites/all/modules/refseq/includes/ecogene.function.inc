<?php 
/*
 * FASTA DNA or Protein format of a sequnce
 *  parameter: 
 *  $sequence the sequence that need to be formatted.
 *  	could be a string or array. 
 *  
 *  	if it is array, will use the main string,
 *  	up string, down string, and the header information.
 *  
 *      if it is a string, it is the main string.
 *  	
 *  $length the character per line. 
 *  $$triplet: if the sequence need have space between triplets. 
 *  			such as: ACT GCA ....
 *  
 */
function fasta_sequence_format($sequence, $length=75, $triplet=FALSE)
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
			 return $FASTA_seq;
	    else 
	    	return $header.$line_terminate.$FASTA_seq;
}
		
			
function calculateMW($sequence,$type='aa')
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


function translationDNAToProtein($sequence)
{
	
	$Codon_Translation_Table = array(

				"GCT" => "A", 
				"GCC" => "A",
				"GCA" => "A",				
				"GCG" => "A",				
				"TGT" => "C",
				"TGC" => "C",				
				"GAT" => "D",
				"GAC" => "D",				
				"GAA" => "E",
				"GAG" => "E",				
				"TTT" => "F",
				"TTC" => "F",				
				"GGT" => "G",
				"GGC" => "G",
				"GGA" => "G",
				"GGG" => "G",								
				"CAT" => "H",
				"CAC" => "H",				
				"ATT" => "I",
				"ATC" => "I",
				"ATA" => "I",				
				"AAA" => "K",
				"AAG" => "K",				
				"TTG" => "L",
				"TTA" => "L",
				"CTT" => "L",
				"CTC" => "L",
				"CTA" => "L",
				"CTG" => "L",				
				"ATG" => "M",				
				"AAT" => "N",
				"AAC" => "N",				
				"CCT" => "P",
				"CCC" => "P",
				"CCA" => "P",
				"CCG" => "P",				
				"CAA" => "Q",
				"CAG" => "Q",				
				"CGT" => "R",
				"CGC" => "R",
				"CGA" => "R",
				"CGG" => "R",
				"AGA" => "R",
				"AGG" => "R",				
				"TCT" => "S",
				"TCC" => "S",
				"TCA" => "S",
				"TCG" => "S",
				"AGT" => "S",
				"AGC" => "S",				
				"ACT" => "T",
				"ACC" => "T",
				"ACA" => "T",
				"ACG" => "T",				
				"GTT" => "V",
				"GTC" => "V",
				"GTA" => "V",
				"GTG" => "V",				
				"TGG" => "W",				
				"NNN" => "X",				
				"TAT" => "Y",
				"TAC" => "Y"
				);
				
		$seqLen = strlen($sequence);
		$protein = "";

		for ($i=0; $i<$seqLen-2;$i+=3)
		{
			$cod = substr($sequence,$i,3);
			if ($cod[0]=='N' || $cod[1]=='N' || $cod[2]=='N')
			{
				$protein = $protein."X";
			}
			else {
				$protein = $protein.$Codon_Translation_Table[$cod];
			}
		}				
		
		return $protein;
}