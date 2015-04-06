<?php 


class RestrictionSite extends DNASequence{
	protected $siteName;
    protected $siteSequence;
  	protected $isPalindrome;
  	protected $siteSequenceRevc;
    

 	public function __construct($sequence,$name=NULL) {
    	if (empty($sequence)) {
      	throw new Exception('Cannot created an instance of RestrictionSite without a specified sequence.');
    	}
    	parent::__construct($sequence);
    	$this->siteName = $name;
    	$this->siteSequence = strtoupper($sequence);
    	$this->setUp();
    	
  	}
	public function siteName() {    	
    	return $this->siteName;	    	
  	}
	public function isPalindrome() {    	
    	return $this->isPalindrome;	    	
  	}
	public function siteSequence() {    	
    	return $this->siteSequence;	    	
  	}
	
  	public function siteSequenceRevc() {    	
    	return $this->$siteSequenceRevc;	    	
  	}
  	
	public function isSeqPalindrome($sequence) {  
		$sequence = strtoupper($sequence);
		$rev_seq = DNASequence:: revc($sequence);
  		
  		$str = str_replace(array('X'), array('N'), $rev_seq);
		$str2 = str_replace(array('X'), array('N'), $sequence);
		  	
    	return ($str == $str2)? TRUE:FALSE;    	
  	}
  	
  	function setUp(){  		
  		
  		$this->isPalindrome = $this->isSeqPalindrome($this->sequence);
  		if($this->isPalindrome){
  			$this->siteSequenceRevc = NULL;
  		}else{
  			$this->siteSequenceRevc =$this->revc($this->siteSequence);
  		}
  		
  	}

	/*
	 * search in the sequence for the site. 
	 * 
	 * $detail is true, return the location 
	 */
	public function search_site($sequence, $detail=FALSE)
	{
		$search = array();
		$replace = array();
		foreach ($this->ambiguity_code() as $key=>$value) {
			$search[] = $key;
			$replace[] = '['.$value.']'
			;
		}
		$matches = array();
		$pattern = "|".str_replace($search , $replace ,$this->siteSequence)."|U";		
		$n_detected = preg_match_all($pattern, strtoupper($sequence), $matches, PREG_OFFSET_CAPTURE | PREG_SET_ORDER);
		
		/*
		 * need to also search the reverse complement if the site is not Palindrome
		 */
		$matches2 = array();
		if(!($this->isPalindrome) && ($this->siteSequenceRevc)){
			$pattern = "|".str_replace($search , $replace ,$this->siteSequenceRevc)."|U";
       		$n_detected += preg_match_all($pattern, strtoupper($sequence), $matches2, PREG_OFFSET_CAPTURE | PREG_SET_ORDER);
       		
       	}
       	if(!$detail){
       		return $n_detected;
       	}else {
       		
       		$detected = array();
       		if($this->isPalindrome){
       			
       			foreach ($matches as $value) {
					$detected[]=array(
						'string' => $value[0][0],
						'offset' => $value[0][1],
						'strand' => '',
					);
				}				
       		}else{
       			foreach ($matches as $value) {
					$detected[]=array(
						'string' => $value[0][0],
						'offset' => $value[0][1],
						'strand' => '+',
					);
				}		
				
       			foreach ($matches2 as $value) {
					$detected[]=array(
						'string' => $value[0][0],
						'offset' => $value[0][1],
						'strand' => '-',
					);
				}		
       		}				
				
			return $detected;
       	}
	}
		
}