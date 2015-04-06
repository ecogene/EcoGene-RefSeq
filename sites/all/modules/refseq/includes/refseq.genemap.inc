<?php 


class GeneMap
{
 
    protected $left, $right;   // the left and the right of the address 
    protected $genome_length; // the length the genome   
    
 	protected $img, $imgfile; // the image object of the map and the file to save the image
 	protected $width, $height;  // width and height of the map
  	protected $margin=array(); //margins for additional inforomation of the map
  	
  	protected $background_color; 
  	protected $fontfile, $color;
	function get_img(){
		return $this->img;
	}
  function get_margin(){
		return $this->margin;
	}

 /**
 * Get details about an image.
 *
 * Drupal supports GIF, JPG and PNG file formats when used with the GD
 * toolkit, and may support others, depending on which toolkits are
 * installed.
 *
 * @param $parametera
 *   keyed array containing information about the gene map:
 *   - "width": Width, of the image area that draws the genes.
 *   - "height": Height, in pixels, of the image area that draws the genes.
 *   - "left": left end of the genome address that will be on the map
 *   - "right": right end of the genome address that will be on the map
 *   	left<right in most of the cases, but the left>right at the address
 *      of begin/end of the genome
 *   - "genome_length": the length of the genome.
 * @param $margin
 *   An array to specify the margin (for title, lengend and etc.) of the image.
 *   0=>left, 0=>right, 0=>top, 0=>bottom, 
 *
 * @return
 *   FALSE, if the file could not be found or is not an image. Otherwise, a
 *   keyed array containing information about the image:
 *   - "width": Width, in pixels.
 *   - "height": Height, in pixels.
 *   - "extension": Commonly used file extension for the image.
 *   - "mime_type": MIME type ('image/jpeg', 'image/gif', 'image/png').
 *   - "file_size": File size in bytes.
 */  
 	public function __construct($parameter=array(), $margin=null) {
    	if (empty($parameter)) {
      	throw new Exception('Cannot created an instance of RestrictionSite without a specified sequence.');
    	}
    	
    	$this->width = $parameter['width'];
    	$this->height = $parameter['height'];
    	$this->left = $parameter['start'];
    	$this->right = $parameter['end'];
    	$this->genome_length = $parameter['genome_length'];
    	
    	$this->margin = is_null($margin)? array('left'=>0,'right'=>0,'top'=>0,'bottom'=>0):$margin;    	
      	$this->fontfile = $parameter['fontfile'];
      	
    	$this->setUp();
  	}
  	function setUp()
  	{
  		$width = $this->width + $this->margin['left'] + $this->margin['right'];
  		$height = $this->height + $this->margin['top'] + $this->margin['bottom'];
  		$this->img = imagecreatetruecolor($width, $height)	or die ("Cannot Initialize new GD image stream");
   		
  		// set the color for different types of genes.
  		imagecolorallocate($this->img,30,144,255);
  		$this->color = array(
  			'Clockwise' => imagecolorallocate($this->img,30,144,255), // clockwise gene
  			'Counterclockwise' => imagecolorallocate($this->img,0, 0, 255), // counter clockwise gene
  			'+' => imagecolorallocate($this->img,30,144,255), // clockwise gene
  			'-' => imagecolorallocate($this->img,0, 0, 255), // counter clockwise gene
  			'rna' => imagecolorallocate($this->img,160,32,240), // rna gene
  			'pseudogene' => imagecolorallocate($this->img,255, 0, 0), // pseudo gene
  			'background' => imagecolorallocate($this->img, 255,255,255), // background color
  			'default' => imagecolorallocate($this->img, 23,43,0),   // default for text and others
  		);
  		if(is_null($this->margin)){
  			
  			$margin = array('left'=>0,'right'=>0,'top'=>0,'bottom'=>0);
  		}
  		imagefill($this->img , 0,0, $this->color['background']);
  		
  	}
 	public function add_lengend($op, $top=0, $left=0, $text_size=8){
 		
 		$text_size = 8;
 		$rect_width = 50;
 		
 		$x = $this->margin['left'] + $left;
 		$y = $this->margin['top'] + $this->height + $top + $text_size;
 		
 		foreach ($op as $name) {
 			
 			foreach ($this->color as $key=>$color) {
 				if(strcasecmp($name, $key)==0){
 					 					
					imagefilledrectangle($this->img, $x, $y, $x+$rect_width, $y+$text_size, $color);
					$x += $rect_width+5;
					
					imagettftext ($this->img, $text_size, 0, $x, $y+$text_size,  $color, $this->fontfile, $name);
					$x +=  strlen($name)*$text_size;
					
 				}
 			}
 		}
 		
 	}
  	public function draw_title($op){
  		switch ($op) {
  			case 'kb': // only the kb
  				$this->draw_title_kb();
  			break;
  			case 'minute': // only the kb
  				$this->draw_title_minute();
  			break;
  			default:
  				$this->draw_title_all();
  				
  			break;
  		}
  		
  	}
/**
 * Draw the kb and/or minute title for the map in the map's top margin area.
 *
 * @param $size
 *   text size.
 * @param $color
 *   An color identifier.
 *
 * @return
 *   FALSE, if the file could not be found or is not an image. Otherwise, a
 *   keyed array containing information about the image:
 *   - "width": Width, in pixels.
 *   - "height": Height, in pixels.
 *   - "extension": Commonly used file extension for the image.
 *   - "mime_type": MIME type ('image/jpeg', 'image/gif', 'image/png').
 *   - "file_size": File size in bytes.
 */  	
  	function draw_title_all($color=null){
  		
  		$bigstep = 1000; // big tic mark at every 1k and the text
  		$smallstep = 200; // small tic mark at every 200 bits
  		
  		$size=7;
  		$line_size = 7;
  		if($this->margin['top']>40){
  			$y = ($this->margin['top'] -40) + 20;
  		}else{
  			$y = 20;
  		}
  		$map_range = ($this->right - $this->left + 1 + $this->genome_length)%$this->genome_length;
  		
  		if(is_null($color)) $color = $this->color['default'];		
		$x_1 = $this->margin['left'] + 1;
		$x_2 = $this->margin['left'] + $this->width;	
		if($this->margin['right']>70){
			$x_2 = $x_2+20;
		}elseif ($this->margin['right']>60){
			$x_2 = $x_2+10;
		}
		
		
		imageline($this->img, $x_1, $y, $x_2, $y, $color);			
		imagettftext ($this->img, $size+1, 0, $x_2+$size, $y-$line_size,  $color, $this->fontfile,'kb' );	
				
		
		$start = ceil($this->left/$bigstep)*$bigstep;
		$stop = ($this->right>$this->left)? $this->right: ($this->right + $this->genome_length);
		
		
		for ($t=$start; $t<=$stop;$t=$t+$bigstep)
		{
			
			$kb_text =  sprintf("%.1f",intval(($t%$this->genome_length))/1000);
			
			$pos_1 = round($this->margin['left']+1 + (($t - $this->left+$this->genome_length)%$this->genome_length)/$map_range*$this->width);
			imageline($this->img, $pos_1, $y-$line_size, $pos_1, $y, $color);	

			imagettftext ($this->img, $size, 0, $pos_1-$line_size, $y-$line_size-2,  $color, $this->fontfile,$kb_text );
			
		}
		
		$start = ceil($this->left/$smallstep)*$smallstep;
		$stop = ($this->right>$this->left)? $this->right: $this->right + $this->genome_length;
			
		for ($t=$start; $t<=$stop;$t=$t+$smallstep)
		{
			$pos_1 = round($this->margin['left'] +1 + (($t - $this->left+$this->genome_length)%$this->genome_length)/$map_range*$this->width);
					
			imageline($this->img, $pos_1,$y-floor($line_size/2),$pos_1, $y,$color);
			
		}
		
		
		// the following draws the minute
		$step_min = 0.0005; // defaut step for minute is 0.5%
		 
		if($this->genome_length*$step_min/$map_range*$this->width<$size*5){
			$step_min = 0.001; 
		}
		$y = $y+1;
		imageline($this->img, $x_1, $y, $x_2, $y, $color);			
		imagettftext ($this->img, $size+1, 0, $x_2+$size, $y+ceil($line_size/2)+$size,  $color, $this->fontfile,'min(Cs)' );	
					
			$bigstep = $step_min*$this->genome_length;
			$smallstep = $bigstep/5;
			
  			$start = ceil($this->left/$bigstep)*$bigstep;
			$stop = ($this->right>$this->left)? $this->right: $this->right + $this->genome_length;
			
	
		
			for ($t=$start; $t<=$stop; $t=$t+$bigstep)
			{
				$min = intval(($t%$this->genome_length))/$this->genome_length*100;
				$min = ($min>100)? ($min-100): $min;
				
				$minute_text =  sprintf("%4.2f",$min);
				$pos_1 = round($this->margin['left']+1 + (($t - $this->left+$this->genome_length)%$this->genome_length)/$map_range*$this->width);
							
				$this->imagelinethick($this->img, $pos_1,1+$y,$pos_1,1+$line_size+$y,$color,1);
				imagettftext ($this->img, $size, 0, $pos_1-$size*2, 5+2*$size+$y,  $color, $this->fontfile,sprintf("%4.2f",$min) );
				
			}
			
			
		
  			$start = ceil($this->left/$smallstep)*$smallstep;
				
			for ($t=$start; $t<=$stop; $t=$t+$smallstep)
			{
				$min = intval(($t%$this->genome_length))/$this->genome_length*100;
				$min = ($min>100)? ($min-100): $min;
				
				$minute_text =  sprintf("%4.2f",$min);
				$pos_1 = round($this->margin['left']+1 + (($t - $this->left+$this->genome_length)%$this->genome_length)/$map_range*$this->width);
							
				$this->imagelinethick($this->img, $pos_1,1+$y,$pos_1,1+ceil($line_size/2)+$y,$color,1);
				
			}
		return ;
  	}
function draw_title_kb($color=null){
  		
  		$bigstep = 1000; // big tic mark at every 1k and the text
  		$smallstep = 200; // small tic mark at every 200 bits
  		
  		$size=7;
  		$line_size = 7;
  		if($this->margin['top']>40){
  			$y = ($this->margin['top'] -40) + 20;
  		}else{
  			$y = 20;
  		}
  		$map_range = ($this->right - $this->left + 1 + $this->genome_length)%$this->genome_length;
  		
  		if(is_null($color)) $color = $this->color['default'];		
		$x_1 = $this->margin['left'] + 1;
		$x_2 = $this->margin['left'] + $this->width;	
		if($this->margin['right']>70){
			$x_2 = $x_2+20;
		}elseif ($this->margin['right']>60){
			$x_2 = $x_2+10;
		}
		
		
		imageline($this->img, $x_1, $y, $x_2, $y, $color);			
		imagettftext ($this->img, $size+1, 0, $x_2+$size, $y-$line_size,  $color, $this->fontfile,'kb' );	
				
		
		$start = ceil($this->left/$bigstep)*$bigstep;
		$stop = ($this->right>$this->left)? $this->right: ($this->right + $this->genome_length);
		
		
		for ($t=$start; $t<=$stop;$t=$t+$bigstep)
		{
			
			$kb_text =  sprintf("%.1f",intval(($t%$this->genome_length))/1000);
			
			$pos_1 = round($this->margin['left']+1 + (($t - $this->left+$this->genome_length)%$this->genome_length)/$map_range*$this->width);
			imageline($this->img, $pos_1, $y-$line_size, $pos_1, $y, $color);	

			imagettftext ($this->img, $size, 0, $pos_1-$line_size, $y-$line_size-2,  $color, $this->fontfile,$kb_text );
			
		}
		
		$start = ceil($this->left/$smallstep)*$smallstep;
		$stop = ($this->right>$this->left)? $this->right: $this->right + $this->genome_length;
			
		for ($t=$start; $t<=$stop;$t=$t+$smallstep)
		{
			$pos_1 = round($this->margin['left'] +1 + (($t - $this->left+$this->genome_length)%$this->genome_length)/$map_range*$this->width);
					
			imageline($this->img, $pos_1,$y-floor($line_size/2),$pos_1, $y,$color);
			
		}
		
	
		return ;
}
  	function draw_title_minute($color=null){
  		
  		$bigstep = 1000; // big tic mark at every 1k and the text
  		$smallstep = 200; // small tic mark at every 200 bits
  		
  		$size=7;
  		$line_size = 7;
  		if($this->margin['top']>20){
  			$y = $this->margin['top'] - 10;
  		}else{
  			$y = 10;
  		}
  		$map_range = ($this->right - $this->left+$this->genome_length + 1)%$this->genome_length;
  		
  		if(is_null($color)) $color = $this->color['default'];		
		$x_1 = $this->margin['left'] + 1;
		$x_2 = $this->margin['left'] + $this->width;	
		
		$y = $size;	
		
  			// the following draws the minute
		$step_min = 0.005; // defaut step for minute is 0.5%
		 
		if($this->genome_length*$step_min/map_range*$this->width<$size*5){
			$step_min = 0.01; 
		}
		$y = $y+1;
		
		imageline($this->img, $x_1, $y, $x_2, $y, $color);			
		imagettftext ($this->img, $size, 0, $x_2+$size*2, $y-ceil($line_size/2),  $color, $this->fontfile,'min(Cs)' );	
		
			
  			$start = ceil($this->left/$step_min)*$step_min;
			$stop = ($this->right>$this->left)? $this->right: $this->right + $this->genome_length;
			
			for ($t=$start; $t<=$stop; $t=$t+$step_min*$this->genome_length)
			{
				$min = intval(($t%$this->genome_length))/$this->genome_length*100;
				$min = ($min>100)? ($min-100): $min;
				
				$minute_text =  sprintf("%4.2f",$min);
				$pos_1 = round($this->margin['left'] + (($t - $this->left)%$this->genome_length)/$map_range*$this->width);
				imageline($this->img, $pos_1, $y-$line_size, $pos_1, $y, $color);	
				imagettftext ($this->img, $size, 0, $pos_1-$line_size, $y,  $color, $this->fontfile,printf("%4.2f",$min) );
				
			}
			$step_min = $step_min/5;
  			$start = ceil($this->left/$step_min)*$step_min;
			$stop = ($this->right>$this->left)? $this->right: $this->right + $this->genome_length;
			
			for ($t=$start; $t<=$stop; $t=$t+$step_min*$this->genome_length)
			{
				$min = intval(($t%$this->genome_length))/$this->genome_length*100;
				$min = ($min>100)? ($min-100): $min;
				
				$minute_text =  sprintf("%4.2f",$min);
				$pos_1 = round($this->margin['left'] + (($t - $this->left)%$this->genome_length)/$map_range*$this->width);
			
				imageline($this->img, $pos_1,$y+floor($line_size/2),$pos_1, $y,$color,1);
			}
		
  	}
  		
		function arrow($im, $x1, $y1, $x2, $y2, $alength, $awidth, $color) {

			$distance = sqrt(pow($x1 - $x2, 2) + pow($y1 - $y2, 2));

			$dx = $x2 + ($x1 - $x2) * $alength / $distance;
			$dy = $y2 + ($y1 - $y2) * $alength / $distance;

			$k = $awidth / $alength;

			$x2o = $x2 - $dx;
			$y2o = $dy - $y2;

			$x3 = $y2o * $k + $dx;
			$y3 = $x2o * $k + $dy;

			$x4 = $dx - $y2o * $k;
			$y4 = $dy - $x2o * $k;
			
			$points = array(
			$x3, $y3,
			$x4, $y4,
			$x2, $y2,

			);
			imagefilledpolygon($im, $points, 3, $color);

		}
		function imagelinethick($image, $x1, $y1, $x2, $y2, $color, $thick = 1)
		{
//			var_dump(array($x1, $y1, $x2, $y2));
			/* this way it works well only for orthogonal lines
			imagesetthickness($image, $thick);
			return imageline($image, $x1, $y1, $x2, $y2, $color);
			*/
			if ($thick == 1) {
				return imageline($image, $x1, $y1, $x2, $y2, $color);
			}
			$t = $thick / 2 - 0.5;
			if($x1 == $x2|| $y1 == $y2)
			{
				
				return imagefilledrectangle($image, round(min($x1, $x2) - $t/3), round(min($y1, $y2) - $t), round(max($x1, $x2) + $t/3), round(max($y1, $y2) + $t), $color);
			}
//			if ($x1 == $x2 || $y1 == $y2) {
//				return imagefilledrectangle($image, round(min($x1, $x2) - $t), round(min($y1, $y2) - $t), round(max($x1, $x2) + $t), round(max($y1, $y2) + $t), $color);
//			}
			
			$k = ($y2 - $y1) / ($x2 - $x1); //y = kx + q
			$a = $t / sqrt(1 + pow($k, 2));
			$points = array(
			round($x1 - (1+$k)*$a), round($y1 + (1-$k)*$a),
			round($x1 - (1-$k)*$a), round($y1 - (1+$k)*$a),
			round($x2 + (1+$k)*$a), round($y2 - (1-$k)*$a),
			round($x2 + (1-$k)*$a), round($y2 + (1+$k)*$a),
			);
			imagefilledpolygon($image, $points, 4, $color);
			return imagepolygon($image, $points, 4, $color);
		}
		
		
/////		
// the elements are draw in the main part of the map
//
//   [$left_margin imx-right_margin] by [top_margin imy-bottom_margin] (width by hight)
//    
//   The elements are supposed to be vertical line and 
//    and their are located from left_end to right_end from the sequence.
//    and the text will be drawn seperated in draw_text function
// 
/////


		function draw_elems($elem = NULL, $para=NULL)
		{
			
			
			$GENOME_LENGTH =  $this->genome_length;
			
			$map_width = $this->width;  //phisycal map horizontal draw range
			$map_range = ($this->right - $this->left+1+$GENOME_LENGTH)%$GENOME_LENGTH; // sequence range in the map
		
			
			$row_space = 40;
			$line_size = 5;
			$text_size = 8;
			
			
			if(isset($para['row_space'])) $row_space = $para['row_space'];
			if(isset($para['line_size'])) $line_size = $para['line_size'];
			if(isset($para['text_size'])) $text_size = $para['text_size'];
			
			$elem_map = array();
			foreach ($elem as $key=>$e) {
				
				if(!empty($e[exon])){
					$line_size = 1;
				}else {
					$line_size = (isset($para['line_size']))? $para['line_size'] : 5;
				}
				$e['start'] = intval($e['start']);
				$e['end'] = intval($e['end']);
				
				
				if(isset($this->color[$e['type']])){
					$color = $this->color[$e['type']];
				}elseif(isset($this->color[$e['strand']])){
					$color = $this->color[$e['strand']];
				}else{
					$color = $this->color['default'];
				}
				
				if(empty($e['row'])) $e['row'] = 1;
				
				$y1 = $this->margin['top']  + ($e['row']-0.5) * $row_space;
				
				if(($e['start'] < $this->left)&& ($e['end'] > $this->left)){
					$e['start'] = $this->left;
				}
				if(($e['end'] > $this->right)&& ($e['start'] < $this->left)){
					$e['end'] = $this->right;
				}
					
				$pos_1 = round($this->margin['left']+1 + (($e['start'] - $this->left + $GENOME_LENGTH)%$GENOME_LENGTH)/$map_range*$map_width);
				$pos_2 = round($this->margin['left']+1 + (($e['end'] - $this->left +$GENOME_LENGTH)%$GENOME_LENGTH)/$map_range*$map_width);
				
				if ($pos_1<($this->margin['left']+1)) $pos_1 = $this->margin['left']+1;
				if ($pos_2>($this->margin['left']+$this->width)) $pos_2 = $this->margin['left']+$this->width;
				
			
				$elem_map[$key] = $e;
				$elem_map[$key]['mappos_left'] = $pos_1;
				$elem_map[$key]['mappos_right'] = $pos_2;
				$elem_map[$key]['mappos_top'] = $y1;
				$elem_map[$key]['mappos_bottom'] = $y1 + $line_size;

				if(!strcmp($e['strand'],'+'))
				{				
	//				
					$this->arrow($this->img, $pos_2-$line_size, $y1, $pos_2, $y1, $line_size, $line_size, $color);
					
					if(($pos_2-$pos_1)>$line_size)
					{
						$this->imagelinethick ($this->img, $pos_1 ,$y1,  $pos_2-$line_size, $y1,$color,$line_size);
					}
					
				}else if(!strcmp($e['strand'],'-'))
				{
				
					
					$this->arrow($this->img, $pos_1+$line_size, $y1, $pos_1, $y1, $line_size, $line_size, $color);
					if(($pos_2-$pos_1)>$line_size)
					{
						$this->imagelinethick ( $this->img,$pos_1+$line_size ,$y1,  $pos_2, $y1,$color,$line_size );
					}
				
				}else if(!strcmp($e['strand'],'+-')) // biodirection
				{
				
					$this->imagelinethick ( $this->img,$pos_1+$line_size ,$y1,  $pos_2-$line_size, $y1,$color,$line_size );
					$this->arrow($this->img, $pos_2-$line_size, $y1, $pos_2, $y1, $line_size, $line_size, $color);
					$this->arrow($this->img, $pos_1+$line_size, $y1, $pos_1, $y1, $line_size, $line_size, $color);
				
				} else {
					
					$this->imagelinethick ( $this->img,$pos_1 ,$y1,  $pos_2, $y1,$color,$line_size );
	
				}
		
				imagettftext($this->img, $text_size, 0, $pos_1, $y1-$line_size-2, $color, $this->fontfile, $e['name']);
				
			}
			return $elem_map;
			
			
		}
		


		function draw_restriction_sites($elem = NULL, $para=NULL)
		{
			
			
			$GENOME_LENGTH =  $this->genome_length;
			
			$map_width = $this->width;  //phisycal map horizontal draw range
			$map_range = ($this->right - $this->left+1+$GENOME_LENGTH)%$GENOME_LENGTH; // sequence range in the map
		
			$row = 1;
			$row_space = 20;
			$line_size = 18;
			$text_size = 8;
			$line_width = 6;
			$color = imagecolorallocate($this->img, rand(0,255),rand(0,255),rand(0,255));
			
			if(isset($para['row_space'])) $row_space = $para['row_space'];
			if(isset($para['line_size'])) $line_size = $para['line_size'];
			if(isset($para['line_width'])) $line_width = $para['line_width'];
			if(isset($para['text_size'])) $text_size = $para['text_size'];
			if(isset($para['color'])) $color = $para['color'];
			if(isset($para['row'])) $row = $para['row'];
			
			$y = $this->margin['top']  + $row * $row_space;
			$this->imagelinethick ($this->img, $this->margin['left']+1 ,$y,  $this->margin['left']+$this->width, $y,$color,1);
			
			$elem_map = array();
			foreach ($elem as $key=>$e) {
				
				$e['start'] = intval($e['start']);
				$e['end'] = intval($e['end']);
			
				$y0 = $this->margin['top'] + 1 + $row* $row_space - $line_size;
				$y1 = $this->margin['top']  + $row * $row_space;
				
				if(($e['start'] < $this->left)&& ($e['end'] > $this->left)){
					$e['start'] = $this->left;
				}
				if(($e['end'] > $this->right)&& ($e['start'] < $this->left)){
					$e['end'] = $this->right;
				}
					
				$pos_1 = round($this->margin['left']+1 + (($e['start'] - $this->left + $GENOME_LENGTH)%$GENOME_LENGTH)/$map_range*$map_width);
//				$pos_2 = round($this->margin['left']+1 + (($e['end'] - $this->left +$GENOME_LENGTH)%$GENOME_LENGTH)/$map_range*$map_width);
				$pos_2 = $line_width + $pos_1;
				
				if ($pos_1<($this->margin['left']+1)) $pos_1 = $this->margin['left']+1;
				if ($pos_2>($this->margin['left']+$this->width)) $pos_2 = $this->margin['left']+$this->width;
				
			
				$elem_map[$key] = $e;
				$elem_map[$key]['mappos_left'] = $pos_1;
				$elem_map[$key]['mappos_right'] = $pos_2;
				$elem_map[$key]['mappos_top'] = $y0;
				$elem_map[$key]['mappos_bottom'] = $y1 ;

				if(!strcmp($e['strand'],'+'))
				{				
	//				
					$this->arrow($this->img, $pos_2-$line_size, ($y1+$y0)/2, $pos_2, ($y1+$y0)/2, $line_size/2, $line_size/2, $color);
			
					
				}else if(!strcmp($e['strand'],'-'))
				{				
					
					$this->arrow($this->img, $pos_1+$line_size, ($y1+$y0)/2, $pos_1, ($y1+$y0)/2, $line_size/2, $line_size/2, $color);
					
				
				} else {
//					var_dump(array($pos_1,$pos_2));
//					$this->imagelinethick ( $this->img,$pos_1 ,$y1,  $pos_2, $y1,$color,$line_size );
					imagefilledrectangle($this->img, $pos_1, $y0, $pos_2, $y1, $color);
	
				}
		
				
				
			}
			return $elem_map;
			
			
		}
		
}



	