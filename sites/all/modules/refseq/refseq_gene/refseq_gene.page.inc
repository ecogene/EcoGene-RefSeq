<?php
function refseq_genepage($accession, $geneid)
{
	global $fna_dir, $faa_dir;
	
	$fnaseqfile=$fna_dir.$accession.'.fna';
	
	if(!is_null(arg(4))){
		$map_magnify = arg(4)>3? 3:arg(4);
		$map_magnify = arg(4)<1? 1:arg(4);
		
		$args = array(
				'genome_acc' => arg(1),
				'geneid' => arg(3),
		 		'map_magnify' => $map_magnify,
				'sites' => explode(',', arg(5)),
				'us' =>(arg(6)=='')?0:intval(arg(6)),
				'ds' =>(arg(7)=='')?0:intval(arg(7)),
				'seq_format'=>(arg(8)=='')?'fasta':arg(8),
		);
	}else{
		$args = array(
			'genome_acc' => arg(1),
			'geneid' => arg(3),
			'map_magnify' => 1,
			'sites' =>array('BamHI','EcoRI','HindIII'),
			'us' =>0,
			'ds' =>0,
			'seq_format'=>'fasta',
		);
	}
	

	// query to get the gene information
	$isCircular = TRUE;
	$range_initial = 10000;
	foreach ($args['sites'] as $key=>$value) {
		$args['sites'][$key] = trim($value)
		;
	}
	$sites_checked =$args['sites'];
	
//	$extra_para = '&'.drupal_http_build_query(drupal_get_query_parameters(array('sites'=>$sites,'map_magnify'=>$map_magnify)) );
		
	
	
	db_set_active('refseq');
	$summary = "summary";
 	$geneTable = "gene";  
 	$refTable = "dbxref";
 	$protTable = "protein";
 	$rnaTable = "rna";
 	$content = "";
 	
 	$GENOME_LENGTH  = db_select($summary, 'summary', array('fetch' => PDO::FETCH_ASSOC))
  						->fields('summary', array('length'))->condition('accession', "$accession%", 'like')
    					->execute()->fetchField();
    $GENOME_INFO  = db_select($summary, 'summary', array('fetch' => PDO::FETCH_ASSOC))
  						->fields('summary')->condition('accession', "$accession%", 'like')
    					->execute()->fetchAssoc();

    $curr_gene = db_select($geneTable, 'gene', array('fetch' => PDO::FETCH_ASSOC))
  			->fields('gene')
    		->condition('accession', $accession)
    		->condition('geneid', $geneid)
     		->execute()->fetchAssoc();
    if(empty($curr_gene))
    	return "No Gene has been found!";

    	
    $query = db_select($refTable, 'dbxref')
  					->condition('accession', $accession)
    					->condition('geneid', $geneid);
    					$query -> condition('type', 'gene');
    			 		$query -> addExpression("CONCAT(dbxref.db, '(', dbxref.ref,')')", 'ref');
    					$ref = $query->execute()->fetchAll(PDO::FETCH_COLUMN);
//
    $curr_gene['ref'] = $ref;
	$curr_gene['name'] = empty($curr_gene['name'])? $curr_gene['locus_tag']:$curr_gene['name'];
	$curr_gene['name'] = empty($curr_gene['name'])? $curr_gene['geneid']:$curr_gene['name'];
	
    
 	$genes = db_select($geneTable, 'gene', array('fetch' => PDO::FETCH_ASSOC))
  			->fields('gene')
    		->condition('accession', $accession)
     		->execute()->fetchAll();

 // GenePage Map information   
    $range = intval($range_initial/pow(2,$args['map_magnify']-1));
    $range_left = $curr_gene['start'] - intval(($range/2));
    $range_right = $curr_gene['start'] + intval($range/2) -1;
   
    if($isCircular){
    	if ($range_left<=0)	$range_left =  $GENOME_LENGTH + $range_left;
    	if ($range_right>$GENOME_LENGTH) $range_right = $range_right % $GENOME_LENGTH;
    }else{
    	if ($range_left<=0)	$range_left =  1;
    	if ($range_right>$GENOME_LENGTH) $range_right = $GENOME_LENGTH;
    }
    
    
    if(!empty($genes)){    
    	$neighbors = array();
    	
     	foreach ($genes as $key => $record){
     		$record['name'] = empty($record['name'])? $record['locus_tag']:$record['name'];
			$record['name'] = empty($record['name'])? $record['geneid']:$record['name'];
	
     		if($range_left<=$range_right){
     			if(   ($record['start']>=$range_left && $record['start']<=$range_right)
     			   || ($record['end']>=$range_left && $record['end']<=$range_right) ){
     				$neighbors[] = $record;
     			}
     			
     		}else{
     			if(  ($record['start']<=$range_right)
     			   ||($record['end']>=$range_left) ){
     				$neighbors[] = $record;
     			}
     			
     		}
     		
    		if($record['geneid']==$geneid){
    				$prev_gene = isset($genes[$key-1])? $genes[$key-1] : $genes[count($genes)-1];  
    				$next_gene = isset($genes[$key+1])? $genes[$key+1] : $genes[0]; 
    				$prev_gene['name'] = empty($prev_gene['name'])? $prev_gene['locus_tag']:$prev_gene['name'];
					$prev_gene['name'] = empty($prev_gene['name'])? $prev_gene['geneid']:$prev_gene['name'];
					$next_gene['name'] = empty($next_gene['name'])? $next_gene['locus_tag']:$next_gene['name'];
					$next_gene['name'] = empty($next_gene['name'])? $next_gene['geneid']:$next_gene['name'];
     		}    		
    	}    	
    }
 	
    // get the intergene inor

   $intergene_info[0] = refseq_intergene_info($prev_gene, $curr_gene,$GENOME_LENGTH);
   $intergene_info[1] = refseq_intergene_info($curr_gene, $next_gene,$GENOME_LENGTH);
   
   //set the genepage title
   drupal_set_title(refseq_genepage_title($curr_gene), PASS_THROUGH);
   
   
	
   $sequence = "";
   DNASequence::refseq_fna_file_read($fnaseqfile, $sequence);

   
   $gene_dna_seq = DNASequence::gene_dna_sequence($sequence,$curr_gene['start'], $curr_gene['end'], $curr_gene['strand'], $args['us'],$args['ds']);
	$gene_dna_seq['header'] = '>GeneID|'.$curr_gene['geneid'].'|ref|'.$GENOME_INFO['accession'].'|'.
   							$curr_gene['start'].'-'.$curr_gene['end'].' '
   							.$GENOME_INFO['taxname'].' '.$GENOME_INFO['replicon'];

   if($curr_gene['type']=='protein' || $curr_gene['type']==''){
   	
   		$faaseqfile=$faa_dir.$accession.'.faa';
   		$gi  = db_select($protTable, 'protein', array('fetch' => PDO::FETCH_ASSOC))
  						->fields('protein', array('protein_id'))
  						->condition('geneid', $curr_gene['geneid'])->condition('Accession', $accession)
    					->execute()->fetchField();
    	if($gi){				
    	
   			$protein_seq = RefSeq::refseq_faa_file_read($faaseqfile, $gi);  
   			$curr_gene['protein_id'] = $gi; 			
    	}   
   }
   // prepare the sequnce information for calculating molecular weight
   if(!empty($protein_seq))	{ $curr_gene['_seq'] = $protein_seq['main_string']; $curr_gene['_type']='aa';}
   elseif($curr_gene['type']=='rna')	{ $curr_gene['_seq'] = $gene_dna_seq['main_string']; $curr_gene['_type']='nt';}
   else { $curr_gene['_seq'] = $gene_dna_seq['main_string']; $curr_gene['_type']='other';}
   $content .= refseq_display_gene_info($curr_gene,$intergene_info,$prev_gene,$next_gene, $GENOME_LENGTH);
   
	
	$range_sequence = ($range_right>$range_left)? substr($sequence, $range_left,$range_right-$range_left+1):substr($sequence, $range_left).substr($sequence, 0, $range_right);
	

    db_set_active("default");
  
    drupal_add_js(drupal_get_path('module', 'refseq_gene') . '/refseq_gene.js');
  	$build['geneinfo'] = array('#markup' => $content,);
  	
  	$build['genemapform'] = drupal_get_form('refseq_genepage_genemap_form', $neighbors,$GENOME_LENGTH, $range_left, $range_right,$args);
//  	$build['sitemap'] = array('#markup' => $sitemap_content,);
  	$build['featuremapform'] = drupal_get_form('refseq_genepage_featuremap_form', $accession,$GENOME_LENGTH, $range_left, $range_right,$args,$neighbors);
  	$build['sitemapform'] =drupal_get_form('refseq_genepage_sitemap_form', $args, $range_sequence,$GENOME_LENGTH, $range_left, $range_right);
  	
    $build['dnaseqform'] = drupal_get_form('refseq_genepage_dnaseq_form', $gene_dna_seq,array_merge(array('start'=>$curr_gene['start'],'end'=>$curr_gene['end'],'strand'=>$curr_gene['strand']), $args));
    if(!empty($protein_seq)) $build['protein'] = drupal_get_form('refseq_genepage_proteinseq_form', $protein_seq);
    
	return $build;
	
}
function refseq_genepage_featuremap_form($form, &$form_state,$accession,$GENOME_LENGTH, $range_left, $range_right,$args,$neighbors){
	 db_set_active("refseq");
	 $feats = db_select('feat', 'feat', array('fetch' => PDO::FETCH_ASSOC))
  			->fields('feat')
    		->condition('accession', $accession)
     		->execute()->fetchAll();
     		
 	if(!empty($feats)){    
    	$features = array();
    	
     	foreach ($feats as $key => $record){
     		
     		$flag_continue = false;
     		
     		foreach ($neighbors as $key_1 => $value_1)
     		{
     			if(($value_1['start']==$record['start']) 
     				&& ($value_1['end']==$record['end'])
     				&& (strcasecmp($value_1['type'], "rna") == 0)
     				)
     					$flag_continue=true;
     		}
     		
     		
     		if((strcasecmp($record['gbkey'], "Src") == 0) || $flag_continue) 
     		{
     			continue;
     		}
     		
     		
     		$record['name'] = empty($record['type'])? $record['type']:$record['gbkey'];
				
     		if($range_left<=$range_right){
     			if(   ($record['start']>=$range_left && $record['start']<=$range_right)
     			   || ($record['end']>=$range_left && $record['end']<=$range_right) ){
     				$features[] = $record;
     			}
     			
     		}else{
     			if(  ($record['start']<=$range_right)
     			   ||($record['end']>=$range_left) ){
     				$features[] = $record;
     			}
     			
     		}
     	}
 	}
     		
	$genemap_content = refseq_gene_featuremap($features, $GENOME_LENGTH, $range_left, $range_right);
	db_set_active("default");
	
	$form['#attributes'] = array('class' => 'genemap-form');
	
	$form['options'] = array(
    '#type' => 'fieldset',
    '#title' => t('FeatureMap'),
    '#collapsible' => TRUE,
    '#collapsed' => !empty($features)? FALSE:TRUE,
    '#weight' => -10,
  
  );
 
  
  $form['options']['genemap'] = array(
    '#markup' => $genemap_content,
  	'#weight' => -30,
  	
  	);
  	  
  	return $form;
}
function refseq_genepage_genemap_form($form, &$form_state,$neighbors,$GENOME_LENGTH, $range_left, $range_right,$args){
	 db_set_active("refseq");
	 foreach ($neighbors as $ind=>$gene_info) {
		 if($gene_info['type']=='protein'){
	   		db_set_active('refseq');
	   		
	   		$result  = db_select('protein', 'protein', array('fetch' => PDO::FETCH_ASSOC))
	  						->fields('protein', array('product','note','start','end','strand'))
	  						->condition('geneid', $gene_info['geneid'])
	    					->execute();
		 	$product = array();
	    	$exon = array();
  
	    	while  ( $row = $result->fetchAssoc()) {
	    		
	    		if(!in_array(($row['product']),($product))){
	    				 $product[] = trim($row['product']);
	    		}
	    		$exon[]= array( 'start' => $row['start'],
	    						'end' => $row['end']);
	    		
	    	}
	    	if($product){
	    			$neighbors[$ind]['product'] = implode("; ", $product);
	    	}
			if(count($exon)>1){
	    			$neighbors[$ind]['exon'] = $exon;
	    	}
	    	db_set_active('default');
   
	 	}
	 	
//	    	if($p){
//	    					$neighbors[$ind]['product'] = $p['product'];
//	    					
//	    	}
//	    	db_set_active('default');
//	   }
		if($gene_info['type']=='rna'){
	   		db_set_active('refseq');
	   		
	   		$result  = db_select('rna', 'rna', array('fetch' => PDO::FETCH_ASSOC))
	  						->fields('rna', array('product','note','start','end','strand'))
	  						->condition('geneid', $gene_info['geneid'])
	    					->execute();
	    	
	    	$product = array();
	    	$exon = array();
  
	    	while  ( $row = $result->fetchAssoc()) {
	    		
	    		if(!in_array(($row['product']),($product))){
	    				 $product[] = trim($row['product']);
	    		}
	    		$exon[]= array( 'start' => $row['start'],
	    						'end' => $row['end']);
	    		
	    	}
	    	if($product){
	    			$neighbors[$ind]['product'] = implode("; ", $product);
	    	}
			if(count($exon)>1){
	    			$neighbors[$ind]['exon'] = $exon;
	    	}
	    	db_set_active('default');
   
	 	}
	 }
	
	$genemap_content = refseq_gene_genemap($neighbors, $GENOME_LENGTH, $range_left, $range_right);
	db_set_active("default");
	
	$form['#attributes'] = array('class' => 'genemap-form');
	
	$form['options'] = array(
    '#type' => 'fieldset',
    '#title' => t('GeneMap'),
    '#collapsible' => TRUE,
    '#collapsed' => FALSE,
    '#weight' => -10,
  
  );
  $form['map_magnify'] = array('#type' => 'value', '#value' => $args['map_magnify']);
  $form['seq_format'] = array('#type' => 'value', '#value' => $args['seq_format']);	
  $form['us'] = array('#type' => 'value', '#value' => $args['us']);
  $form['ds'] = array('#type' => 'value', '#value' => $args['ds']);
  $form['sites'] = array('#type' => 'value', '#value' => $args['sites']);
  
  $form['options']['genemap'] = array(
    '#markup' => $genemap_content,
  	'#weight' => -30,
  	
  	);
  	 $form['options']['maganify_plus'] = array(
    	'#type' => 'image_button',
  		'#src' => drupal_get_path('module', 'refseq').'/files/mag_plus.png',
		'#submit' => array('zoom_in_button_callback','refseq_genepage_genemap_form_submit'),
  	  	'#attributes' => array(
    		'class' => array('magnify'), // change to just 'text' for Drupal 6
  	  	),
  	
  	); 	
  $form['options']['maganify_min'] = array(
    '#type' => 'image_button',
  	'#submit' => array('zoom_out_button_callback','refseq_genepage_genemap_form_submit'),
  	'#src' => drupal_get_path('module', 'refseq').'/files/mag_min.png',
  	'#attributes' => array(
    		'class' => array('magnify'), // change to just 'text' for Drupal 6
  	  	),
  	); 	
  
  	return $form;
}
function zoom_in_button_callback($form, &$form_state){
	
	$map_magnify = intval($form_state['values']['map_magnify']);

	$map_magnify++;
	$form_state['values']['map_magnify'] = $map_magnify>3? 3:$map_magnify;	
	
	return $form;
}
function zoom_out_button_callback($form, &$form_state){
	$map_magnify = intval($form_state['values']['map_magnify']);	 
	$map_magnify--;
	$form_state['values']['map_magnify'] = $map_magnify<1? 1:$map_magnify;		

	return $form;
}


function refseq_genepage_genemap_form_submit($form, &$form_state){

	$map_magnify = intval($form_state['values']['map_magnify']);

	 $args = array(
    	$map_magnify,
    	implode(',', array_filter($form_state['values']['sites'])),
    	$form_state['values']['us'],
    	$form_state['values']['ds'],    
    	$form_state['values']['seq_format'],      	
    	);

  $form_state['redirect'] = 'refseq/'.arg(1).'/gene/' .arg(3). '/'. implode('/', $args);
}
function refseq_genepage_sitemap_form($form, &$form_state,$args, $range_sequence,$GENOME_LENGTH, $range_left, $range_right){
	

	 db_set_active("refseq");
	$sitemap_content = refseq_gene_sitemap( $args['sites'], $range_sequence,$GENOME_LENGTH, $range_left, $range_right);
	 db_set_active("default");
	$form['#attached'] = array(
    	'css' => array   (
        	'type' => 'file',
        	'data' => drupal_get_path('module', 'refseq_gene') . '/refseq_gene.css',
    	),
    
	);
  $form['#attributes'] = array('class' => 'genemap-form');	
  $form['map_magnify'] = array('#type' => 'value', '#value' => $args['map_magnify']);
  $form['seq_format'] = array('#type' => 'value', '#value' => $args['seq_format']);	
  $form['us'] = array('#type' => 'value', '#value' => $args['us']);
  $form['ds'] = array('#type' => 'value', '#value' => $args['ds']);
	
  $form['options'] = array(
    '#type' => 'fieldset',
    '#title' => t('Restriction sites map'),
    '#collapsible' => TRUE,
    '#collapsed' => FALSE,
    '#weight' => -10,
  
  );
  $form['options']['para'] = array(
     '#type' => 'container',
	 '#theme' => 'refseq_gene_inline_form',
  	  '#attributes' => array(
    	'class' => array('refseq_genepage_inline_form'), // change to just 'text' for Drupal 6
  	  ),
  
	  '#weight' => -8,
  );
 $form['options']['para']['sites'] = array(
  		'#type' => 'textfield',
  		'#default_value' => implode(',', $args['sites']),
// 		'#description' => 'Input Enzymes here.',
  		'#size' => 80,
		'#field_prefix' =>l('Pick Additional Restriction Enzymes: ', 'refseq/sitesinfo',array('attributes' =>  array('title'=>t('Full list of enzymes')))),
  		'#weight' => -8,
	);
	

  
	$form['options']['para']['actions'] = array('#type' => 'actions', '#weight' => 0,);
 	$form['options']['para']['actions']['submit'] = array(
    '#type' => 'submit',
    '#value' => t('Submit'),
  	'#weight' => 0,
  	);
  
  	$form['options']['sitemap'] = array(
    '#markup' => '<p>'.$sitemap_content.'</p>',
  	'#weight' => -30,
  	
  	);
  	
  	return $form;
  	

}

function refseq_genepage_sitemap_form_submit($form, &$form_state)
{
	
	 $args = array(
    	$form_state['values']['map_magnify'],
    	$form_state['values']['sites'],
//    	implode(',', array_filter($form_state['values']['sites'])),
    	$form_state['values']['us'],
    	$form_state['values']['ds'],    
    	$form_state['values']['seq_format'],      	
    	);

  $form_state['redirect'] = 'refseq/'.arg(1).'/gene/' .arg(3). '/'. implode('/', $args);
}

function refseq_genepage_proteinseq_form($form, &$form_state,$protein_seq)
{
	$form['#attached'] = array(
    	'css' => array   (
        	'type' => 'file',
        	'data' => drupal_get_path('module', 'refseq_gene') . '/refseq_gene.css',
    	),
    
	);
	
	$string =  DNASequence::fasta_sequence_format($protein_seq);	
	$form['options'] = array(
    '#type' => 'fieldset',
    '#title' => t('Protein sequence'),
    '#collapsible' => TRUE,
    '#collapsed' => FALSE,
    '#weight' => -10,
  
  	);
  
  $form['options']['sequence'] = array(
    '#markup' => "<pre>$string</pre>",
  	'#weight' => 30,
    
  	);
  	return $form;
}
function refseq_genepage_dnaseq_form_submit($form, &$form_state)
{
	
	 $args = array(
    	$form_state['values']['map_magnify'],
    	implode(',', array_filter($form_state['values']['sites'])),
    	$form_state['values']['us'],
    	$form_state['values']['ds'],    
    	$form_state['values']['seq_format'],      	
    	);

  $form_state['redirect'] = 'refseq/'.arg(1).'/gene/' .arg(3). '/'. implode('/', $args);
}
function theme_refseq_gene_inline_form($variables) {
  $form = $variables['form'];

  $output = '<table><tr>';
  foreach (element_children($form) as $element) {
    $output .= '<td width="150px">' . drupal_render($form[$element]) . '</td>';
  }
  $output .= '</tr></table>';

  return $output;
}
function refseq_genepage_dnaseq_form($form, &$form_state,$gene_dna_seq,$args)
{
	
	$us = $args['us'];
	$ds = $args['ds'];
	
	if(empty($us) && empty($ds)){		
		$fasta_string = DNASequence::fasta_sequence_format($gene_dna_seq);			
	}else{
		$fasta_string = DNASequence::fasta_sequence_format($gene_dna_seq, 57,TRUE);
	}
		$coor_string = DNASequence::coor_sequence_format($gene_dna_seq, $args);
	
	$form['#attached'] = array(
    	'css' => array   (
        	'type' => 'file',
        	'data' => drupal_get_path('module', 'refseq_gene') . '/refseq_gene.css',
    	),
    
	);
	
  $form['map_magnify'] = array('#type' => 'value', '#value' => $args['map_magnify']);
  $form['sites'] = array('#type' => 'value', '#value' => $args['sites']);
	
  $form['options'] = array(
    '#type' => 'fieldset',
    '#title' => t('DNA sequence'),
    '#collapsible' => TRUE,
    '#collapsed' => FALSE,
    '#weight' => -10,
  
  );
  $form['options']['para'] = array(
     '#type' => 'container',
	 '#theme' => 'refseq_gene_inline_form',
  	  '#attributes' => array(
    	'class' => array('refseq_genepage_inline_form'), // change to just 'text' for Drupal 6
  	  ),
  
	  '#weight' => -8,
  );
 $form['options']['para']['us'] = array(
  		'#type' => 'textfield',
  		'#default_value' => $us,
  		'#size' => 4,
		'#field_prefix' => t('Upstream: '),
  		'#field_suffix' => t(' bp'),
  		'#weight' => -8,
	);
	$form['options']['para']['ds'] = array(
  		'#type' => 'textfield',
  		'#default_value' => $ds,
  		'#size' => 4,
		'#field_prefix' => t('Downstream: '),
  		'#field_suffix' => t(' bp'),
  		'#weight' => -8,
	);
	
 

  
	$form['options']['para']['actions'] = array('#type' => 'actions', '#weight' => 0,);
 	$form['options']['para']['actions']['submit'] = array(
    '#type' => 'submit',
    '#value' => t('Submit'),
  	'#weight' => 0,
  	);
  	$form['options']['para']['seq_format'] = array(
    '#type' => 'radios',
    '#options' => array(
 					'fasta'=> 'FASTA',
 					'coor'=> 'Coordinates',),
    '#default_value' => $args['seq_format'],
    
  );
  	$form['options']['fasta'] = array(
  	 '#type' => 'container',
  	'#states' => array(
      'visible' => array(
        'input[name="seq_format"]' => array('value' => 'fasta'),
      ),
      ),
  	);
  	$form['options']['fasta']['sequence'] = array(
    '#markup' => "<pre>$fasta_string</pre>",
  	'#weight' => 30,
  	
  	);
  	
  	$form['options']['coor'] = array(
  	 '#type' => 'container',
  	'#states' => array(
      'visible' => array(
        'input[name="seq_format"]' => array('value' => 'coor'),
      ),
      ),
  	);
  	$form['options']['coor']['sequence'] = array(
    '#markup' => "<pre>$coor_string[0]</pre><pre class='sequence'>$coor_string[1]</pre><pre class='coorrel'>$coor_string[2]</pre>",
  	'#weight' => 30,
  	
  	);
  
  	return $form;
}
function refseq_genepage_title($gene_info)
{
	$name = empty($gene_info['name'])? $gene_info['locus_tag']:$gene_info['name'];
	$name = empty($name)? $gene_info['geneid']:$name;
	$title= '<h1 id="gene_page_title"> GenePage for the <span id="g_name">'.$name;
	
	if($gene_info['pseudo']=="true")
		$title .=	'\'</span><span id="pseudo_tag"> pseudogene </span> ';	
	else
		$title .= '</span> gene ';
//	$title .=  'of <span id="escherichia">Escherichia coli</span> K-12 </h1>';
	$title .=  '</h1>';
	return $title;
}
function refseq_intergene_info($row_left_gene, $row_right_gene, $GENOME_LENGTH)
{
	$intgene_name = $row_left_gene['name'].'_'.$row_right_gene['name'];
	
	if(($row_left_gene['strand']=='+') && ($row_right_gene['strand']=='+'))
	{
		$intgene_orientation = 'Codirectional+';
	}
	elseif (($row_left_gene['strand']=='+') && ($row_right_gene['strand']=='-'))
	{
		$intgene_orientation = 'Convergent';
	}
	elseif (($row_left_gene['strand']=='-' )&& ($row_right_gene['strand']=='+'))
	{
		$intgene_orientation = 'Divergent';
	}
	else
	{
		$intgene_orientation = 'Codirectional-';
	}
	
	if($row_right_gene['start']>$row_left_gene['end']+1){
		$overLap = 0;
		$left_end = $row_left_gene['end'] + 1;
		$right_end = $row_right_gene['start'] - 1;
		$intgene_length = $right_end - $left_end + 1;
		$intgene_length = $intgene_length.' bp gap';
		$centisome = $left_end / $GENOME_LENGTH * 100;

	}
	elseif($row_right_gene['start']<$row_left_gene['end']+1){
		$overLap = 1;
		$left_end = $row_right_gene['start'];
		$right_end = $row_left_gene['end'];
		$intgene_length = $right_end - $left_end + 1;
		if($intgene_length>($GENOME_LENGTH/2))//beginning of the genome
		{
			$left_end = $row_left_gene['end']+1;
			$right_end = $row_right_gene['start']-1;
			$intgene_length = $right_end - $left_end + 1 + $GENOME_LENGTH;
			$intgene_length = $intgene_length.' bp gap';
		}else 
		{
			$intgene_length = $intgene_length.' bp overlap';
		}	
		
		$centisome = $left_end / $GENOME_LENGTH * 100;

	}else {
		$intgene_length = '0 bp';
		$left_end = 'Null';
		$right_end = 'Null';
		$centisome = $row_right_gene['start'] / $GENOME_LENGTH * 100;
		
	}
	$centisome = sprintf("%.2f",$centisome);
	$intergene_info['name'] = $intgene_name;
	$intergene_info['length'] = $intgene_length;
	$intergene_info['ori']= $intgene_orientation;
	$intergene_info['left_end'] = $left_end;
	$intergene_info['right_end'] = $right_end;
	$intergene_info['cent'] = $centisome;
	return $intergene_info;
}



function refseq_display_gene_info($gene_info,$intergene_info,$prev_gene,$next_gene, $GENOME_LENGTH)
{
	global $base_url;
	
	$sites = isset($_GET['sites'])? $_GET['sites']:'BamHI, EcoRI, HindIII';

	$sites_array = explode(', ', $sites);
	$map_magnify = isset($_REQUEST['map_magnify'])? $_REQUEST['map_magnify']:1;	
	$extra_para = '&'.drupal_http_build_query(drupal_get_query_parameters(array('sites'=>$sites,'map_magnify'=>$map_magnify)) );
	
 // this part of the code deals the join statement for protein and rna	
	if($gene_info['type']=='protein'){
	   		db_set_active('refseq');
	   		
	   		$result  = db_select('protein', 'protein', array('fetch' => PDO::FETCH_ASSOC))
	  						->fields('protein', array('product','note','start','end','strand'))
	  						->condition('geneid', $gene_info['geneid'])
	    					->execute();
		 	$product = array();
	    	$note = array();  
	    	$exon = array();
	    	while  ( $row = $result->fetchAssoc()) {
	    			    		
	    		if(!in_array(($row['product']),($product))){
	    				 $product[] = trim($row['product']);	    				 
	    		}	
	    		if(!in_array(($row['note']),($note))){
	    				
	    				 $note[] = trim($row['note']);
	    		}	
	    		$exon[]= array( 'start' => $row['start'],
	    						'end' => $row['end']);    		
	    	}
	    	if(!empty($product)){
	    			$product_name = implode("; ", $product);
	    	}
			if(!empty($note)){
	    			$product_note = implode("; ", $note);
	    	}
			if(count($exon)>1){
	    			$gene_info['exon'] = $exon;
	    	}
	    	db_set_active('default');
   
	 }

	if($gene_info['type']=='rna'){
	   		db_set_active('refseq');
	   		
	   		$result  = db_select('rna', 'rna', array('fetch' => PDO::FETCH_ASSOC))
	  						->fields('rna', array('product','note','start','end','strand'))
	  						->condition('geneid', $gene_info['geneid'])
	    					->execute();
	    	
			$product = array();
	    	$note = array();  
			$exon = array();
  
	       	while  ( $row = $result->fetchAssoc()) {
	    			    		
	    		if(!in_array(($row['product']),($product))){
	    				 $product[] = trim($row['product']);	    				 
	    		}	
	    		if(!in_array(($row['note']),($note))){
	    				
	    				 $note[] = trim($row['note']);
	    		}	 
	    		$exon[]= array( 'start' => $row['start'],
	    						'end' => $row['end']);   		
	    	}
	    	if(!empty($product)){
	    			$product_name = implode("; ", $product);
	    	}
			if(!empty($note)){
	    			$product_note = implode("; ", $note);
	    	}
			if(count($exon)>1){
	    			$gene_info['exon'] = $exon;
	    	}
	    	db_set_active('default');   
	 }
	$join_text = "";	
	if(!empty($gene_info['exon'])){
			
			$address_text = array();
			foreach ($gene_info['exon'] as $value) {
				$address_text[]=$value['start'].",".$value['end'];
			}
			$address_text = implode("...", $address_text);
			if($gene_info['strand']=='-'){
				$address_text = "complement(join(".$address_text."))";
			}else{
				$address_text = "join(".$address_text.")";
			}
	}
	
	
	if(!empty($address_text)){
	 	$join_text = '<span title="'.$address_text.'" style ="border-bottom: thin dotted; color: #0070B3;"> join</span>';
	}

	$name = empty($gene_info['name'])? $gene_info['locus_tag']:$gene_info['name'];
	$content = '<div id="gene_info">';
	
	$content .= '
	<table>
		<tr>
			<td>
				Primary Gene Name: <strong>'.$name;
	$content.=	(strcmp($gene_info['pseudo'],"true")==0)?'\'':'';
	$content.=	'</strong><br /> 
				GeneID: '. $gene_info['geneid'].
				'<br />Gene Page Link: 
				<a href="http://www.ncbi.nlm.nih.gov/gene/'.$gene_info['geneid'].'"> NCBI </a>';
				if(!empty($gene_info['protein_id'])){
					$content.='<a href="http://www.uniprot.org/uniprot/?query='.$gene_info['protein_id'].'"> UniProtKB </a>';
				}
				$content.='<br />  
				Locus Tag: '. $gene_info['locus_tag'].'<br />';
				if(!empty($gene_info['ref']))
					$content.= 'Other Gene Identifier: '.implode('; ', $gene_info['ref']);
	
	if(!empty($product_name) || !empty($product_note)){
		$content.=	'<br />';
		if(!empty($product_name)) 	
			$content.=	'<br /> Gene Product: '.$product_name;
		if(!empty($product_note)) 	
			$content.=	'<br /> Note: '.$product_note;
		
	}
	$content.=	'</td>
		</tr></table>	';
	
	
				
	$content .= 
	'<table>
		<tr>
			<td>
				<A HREF="'.$base_url.'?q=refseq/'.$gene_info['accession'].'/gene/'.
				$prev_gene['geneid'].$extra_para.'"><IMG SRC="'.$base_url.'/sites/all/images/cylarrw.gif" 
				title="Previous gene by address"></A>
			</td>
			<td>
				Left End: '.$gene_info['start'].' <br />
				<span class="intergene" onmouseover="toggleDiv(\'left_intergene\', 1);" onmouseout="toggleDiv(\'left_intergene\', 0);")">Left Intergenic Region</span>
				<div id="left_intergene" class="intergene">
					<p> Name: '.$intergene_info[0]['name'].'</p>
					<p> Length: '.$intergene_info[0]['length'].'</p>
					<p> Orientation: '.$intergene_info[0]['ori'].'</p>
					<p> Left_end: '.$intergene_info[0]['left_end'].'</p>
					<p> Right_end: '.$intergene_info[0]['right_end'].'</p>
					<p> Centisome: '.$intergene_info[0]['cent'].'</p>
				</div>
			</td>
			<td>
				Genomic Address '.$join_text.'<br />'.
				($gene_info['strand']=='+'? 'Clockwise' :'Counterclockwise') .'<br />'.
				'Minute or Centisome (%) = '.number_format($gene_info["start"] / $GENOME_LENGTH * 100,2) .'<br />'.
				(($gene_info["_type"]=='aa'||$gene_info["_type"]=='rna')?
				'Molecular Weight: '.number_format(DNASequence::calculateMW($gene_info["_seq"],$gene_info["_type"]),2).
				'------'.DNASequence::calculate_length($gene_info["_seq"],$gene_info["_type"]):'').'<br />'.
			'</td>
			<td>
				Right End: '.$gene_info['end'].' <br />
				<span class="intergene" onmouseover="toggleDiv(\'right_intergene\', 1);" onmouseout="toggleDiv(\'right_intergene\', 0);")">Right Intergenic Region</span>
				<div id="right_intergene" class="intergene">
				<p> Name: '.$intergene_info[1]['name'].'</p>
					<p> Length: '.$intergene_info[1]['length'].'</p>
					<p> Orientation: '.$intergene_info[1]['ori'].'</p>
					<p> Left_end: '.$intergene_info[1]['left_end'].'</p>
					<p> Right_end: '.$intergene_info[1]['right_end'].'</p>
					<p> Centisome: '.$intergene_info[1]['cent'].'</p>
				</div>
			</td>	
			<td>
				<A HREF="'.$base_url.'?q=refseq/'.$gene_info['accession'].'/gene/'.
				$next_gene['geneid'].$extra_para.'"><IMG SRC="'.$base_url.'/sites/all/images/cyrarrw.gif" 
				title="Next gene by address"></A>
				
			</td>
		</tr>
	</table>';
	$content .= '</div>';
						

	return $content;
}
function refseq_gene_genepage_form($form, &$form_state, $values, $statuses=FALSE)
{
	
}
function refseq_gene_sitemap($sites, $sequence,$GENOME_LENGTH, $left_end, $right_end){
	$content = "";
//	$tmpfname="temp/"."sites".rand().".png";
	global $base_url;
	
	$margin = array('left'=>20,'right'=>80,'top'=>10,'bottom'=>10);
	$fontfile = drupal_get_path('module', 'refseq'). "/files/verdana.ttf"; 
	
	$row_space=18;
    $img_width = 850;
    $text_size = 8;
    $line_size= 15;
    db_set_active('refseq'); 

	$enzyme_tbl = "t_enzyme";
	
	$query = db_select($enzyme_tbl,'enzyme')->fields('enzyme');
	$query->condition('name',$sites,'IN');
	$records = $query->orderBy('name')->execute()->fetchAll();;
	
	$sites_obj = array();
	
	foreach ($records as $value) {
	
		$enzyme_name =trim($value->name);
		$pattern_disply = trim($value->sequence);
		$is_symmetrical = $value->is_symmetrical;
		$pattern = trim($value->seq);
	
		
		$sites_obj[] = new RestrictionSite($pattern,$enzyme_name);
		
	}
	
	 $img_high =  count($sites_obj)*$row_space; 	 
  
	$para = array(	
		'width'=>$img_width,
    	'height'=>$img_high,
    	'start'=>$left_end,
    	'end'=>$right_end,
    	'genome_length'=>$GENOME_LENGTH,
    	'fontfile'=>$fontfile,
	);
     
	$site_map = new GeneMap($para,$margin);
	
	$img = $site_map->get_img();
	$margin = $site_map->get_margin();
	$row = 1;
	$map_site_map = array();
	$site_map_para= array(
			'row_space' => $row_space,
			'line_size' => $line_size,
			'text_size' => $text_size,
			'line_width' => $line_size/6,
			);
	
	foreach ($sites_obj as $obj) {
		$detected = $obj->search_site($sequence, TRUE); // search sites for given sequence
		
		foreach ($detected as $key=>$value) {
			$detected[$key]['start'] = ($value['offset']+$left_end)%$GENOME_LENGTH;
			$detected[$key]['end'] = ($value['offset']+$left_end)%$GENOME_LENGTH;
			
		}
		
		$site_map_para ['row'] = $row;
		$site_map_para ['color'] =   imagecolorallocate($img,rand(0,255),rand(0,255),rand(0,255));
		$s_map = $site_map->draw_restriction_sites($detected,$site_map_para);
		
		foreach ($s_map as $value) {
			$value['start'] = ($value['offset']+1+$left_end)%$GENOME_LENGTH;
			$value['end'] = ($value['offset']+1+$left_end)%$GENOME_LENGTH;
			$map_site_map[] = $value;
		}
		$pos_x = $margin['left']+ $img_width+2;
		$pos_y = $margin['top']  + $row * $row_space;

		imagettftext($img, $text_size, 0, $pos_x, $pos_y, $site_map_para ['color'], $fontfile, $obj->siteName());
		$map_site_map[] = array( 'string' => 'reference from REBASE',
								'url' => 'http://rebase.neb.com/rebase/enz/'.$obj->siteName().'.html',
								'mappos_left' => $pos_x,
      							'mappos_right' => $pos_x + strlen($obj->siteName())*$text_size,
      							'mappos_top' => $pos_y -$text_size,
      							'mappos_bottom' => $pos_y,
								'start'=>'');
		$row++;
	}
	
	
	$tmpfname = drupal_tempnam('public://', 'site').'.png';
		
	Imagepng($img,drupal_realpath($tmpfname));
	
	ImageDestroy($img);

	$content.= '<map name="area_site_map">';
	$link = $base_url.'?q=refseq/'.arg(1).'/gene/';
	
	foreach ($map_site_map as $t=>$value) {
		$link = isset($value['url'])? $value['url']: 'javacsript:';
		$content.= "<AREA NAME=\"area_site_map".$t."\" 
				COORDS=\"".$value['mappos_left'].",".$value['mappos_top'].",".$value['mappos_right'].",".$value['mappos_bottom']."\" 
				title=\"".$value['start']." ".$value['string']."\" target=\"_parent\" HREF=\"".$link."\"".">";
	}
		
		$content.= "</map>";
//		$content.= '<div class="EG_CollapsibleMap_TitleDiv">
//					<div>
//						<a class="EG_CollapsibleMap_TitleAhref" href="javascript:void(0)" onclick="toggleBlockDiv(\'sitemap\', \'sites_CollapsibleMap\');">
//						<img id="sites_CollapsibleMap" class ="cl_CollapsibleArea_collapsing" src="sites/all/images/empty.gif" border="1">
//						<span class="EG_CollapsibleMap_Title">SitesMap</span></a>
//						<div class="EG_CollapsibleMap_HrDiv"> <hr class="EG_CollapsibleMap_Hr"></div>
//					</div></div>';
		
		$content.= '<div id="sitemap">';	
		$content.= '<img style="display:block;" SRC="'.file_create_url($tmpfname).'" usemap="#area_site_map" />';
		$content.='</div>';
		
	return $content;	
	
}
function refseq_gene_featuremap($genes, $GENOME_LENGTH, $left_end, $right_end){
	 
	

	 global $base_url;
	 $content = '';
	$row_space=40;
    $img_width = 850;
    $img_high =  2*$row_space; 
    $fontfile = drupal_get_path('module', 'refseq'). "/files/verdana.ttf"; 
    
     $text_size = 8;
     $line_size = 5;
     
	$para = array(	
		'width'=>$img_width,
    	'height'=>$img_high,
    	'start'=>$left_end,
    	'end'=>$right_end,
    	'genome_length'=>$GENOME_LENGTH,
    	'fontfile'=>$fontfile,
	);
     $margin = array('left'=>20,'right'=>80,'top'=>40,'bottom'=>30);
       	
	$gene_map = new GeneMap($para,$margin);

	if($GENOME_LENGTH<1000000)	{
		$gene_map->draw_title('kb');
	}
	else {
		$gene_map->draw_title();
	}
	
	foreach ($genes as $key=>$e) {
				$genes[$key]['row'] = ($key+1)%2+1;
				if($e['pseudo']=='true')	$genes[$key]['type'] = 'pseudogene';
	}
	
	$mapofsites = array();
	$mapofsites = $gene_map->draw_elems($genes);
	
//	$gene_map->add_lengend(array('misc_RNA','nc_RNA','rRNA','tmRNA'));
//	var_dump($para);
//	var_dump($genes);

	$tmpfname = drupal_tempnam('public://', 'feature').'.png';
		
	$img = $gene_map->get_img();
	Imagepng($img,drupal_realpath($tmpfname));
	ImageDestroy($img);

	$content.= '<map name="mapfeaturemap">';
//	var_dump($mapofsites);
	foreach ($mapofsites as $t=>$value) {
		$content.= "<AREA NAME=\"area_feature".$t."\" 
				COORDS=\"".$value['mappos_left'].",".$value['mappos_top'].",".$value['mappos_right'].",".$value['mappos_bottom']."\" 
				title=\"".$value['start']." - ".$value['end']."\n".$value['note']."\" target=\"_parent\" HREF=\"".$link."\"".">";
	}
		
		$content.= "</map>";

		$content .= '<img style="display:block;" SRC="'.file_create_url($tmpfname).'" usemap="#mapfeaturemap" />';
	return $content;	
	
}
function refseq_gene_genemap($genes, $GENOME_LENGTH, $left_end, $right_end){
	 
	
	 
//	$tmpfname="temp/".$genes['geneid']."sites".rand().".png";
	 global $base_url;
	 $content = '';
	$row_space=40;
    $img_width = 850;
    $img_high =  2*$row_space; 
    $fontfile = drupal_get_path('module', 'refseq'). "/files/verdana.ttf"; 
    
     $text_size = 8;
     $line_size = 5;
     
	$para = array(	
		'width'=>$img_width,
    	'height'=>$img_high,
    	'start'=>$left_end,
    	'end'=>$right_end,
    	'genome_length'=>$GENOME_LENGTH,
    	'fontfile'=>$fontfile,
	);
     $margin = array('left'=>20,'right'=>80,'top'=>40,'bottom'=>30);
       	
	$gene_map = new GeneMap($para,$margin);

	if($GENOME_LENGTH<1000000)	{
		$gene_map->draw_title('kb');
	}
	else {
		$gene_map->draw_title();
	}
	$genes_display = $genes;
	foreach ($genes as $key=>$e) {
				$genes[$key]['row'] = ($key+1)%2+1;
				if($e['pseudo']=='true')	$genes[$key]['type'] = 'pseudogene';
	}
	$ind = 0;
	foreach ($genes as $key=>$e) {
		
		if(!empty($genes[$key]['exon'])){
			
			$address_text = array();
			foreach ($genes[$key]['exon'] as $value) {
				$address_text[]=$value['start'].",".$value['end'];
			}
			$address_text = implode("...", $address_text);
			if($genes[$key]['strand']=='-'){
				$address_text = "complement(join(".$address_text."))";
			}else{
				$address_text = "join(".$address_text.")";
			}
			$genes_display[$ind] = $genes[$key];
			
			$genes_display[$ind]['address_text'] = $address_text;
			$genes_display[$ind]['row'] = ($key+1)%2+1;
			if($genes[$key]['pseudo']=='true'){
				$genes_display[$ind]['type'] = 'pseudogene';
			}	
			$genes_t = $genes_display[$ind];
			
			$genes_display[$ind]['strand']='';
			$genes_display[$ind]['name']='';
			$ind++;
			$is_name_set = false;
			foreach ($genes_t['exon'] as $value) {
				
				$genes_display[$ind] = $genes_t;
				unset($genes_display[$ind]['exon']);
				$genes_display[$ind]['start'] = $value['start'];
				$genes_display[$ind]['end'] = $value['end'];
				if($is_name_set)
					$genes_display[$ind]['name'] = '';
				$is_name_set = true;
				$ind++;
			}
		}else{
				$genes_display[$ind] = $genes[$key];
				$genes_display[$ind]['row'] = ($key+1)%2+1;
				if($e['pseudo']=='true'){
					$genes_display[$ind]['type'] = 'pseudogene';
				}	
				$genes_display[$ind]['address_text'] = $genes[$key]['start']."...".$genes[$key]['end'];
				$ind++;
		}
	}
	
	$mapofsites = array();
//	$mapofsites = $gene_map->draw_elems($genes);
	$mapofsites = $gene_map->draw_elems($genes_display);
//	var_dump($genes);
//	var_dump($genes_display);
	$gene_map->add_lengend(array('RNA','Pseudogene','Clockwise','Counterclockwise'));
	
//	$tmpfname="temp/gene".$genes['geneid'].rand().".png";
	$tmpfname = drupal_tempnam('public://', 'gene').'.png';
		
	$img = $gene_map->get_img();
	Imagepng($img,drupal_realpath($tmpfname));
	ImageDestroy($img);
	
	$content.= '<map name="mapgenemap">';
	$link = $base_url.'?q=refseq/'.arg(1).'/gene/';
	foreach ($mapofsites as $t=>$value) {
		$content.= "<AREA NAME=\"area_gene".$t."\" 
				COORDS=\"".$value['mappos_left'].",".$value['mappos_top'].",".$value['mappos_right'].",".$value['mappos_bottom']."\" 
				title=\"Product: ".$value['product']."&#13;Address: ".$value['address_text']."\" target=\"_parent\" HREF=\"".$link.$value['geneid']."\"".">";
		
//		$content.= "<AREA NAME=\"area_gene".$t."\" 
//				COORDS=\"".$value['mappos_left'].",".$value['mappos_top'].",".$value['mappos_right'].",".$value['mappos_bottom']."\" 
//				title=\"Product: ".$value['product']."&#13;Address: ".$value['start']."... ".$value['end']."\" target=\"_parent\" HREF=\"".$link.$value['geneid']."\"".">";
	}
		
		$content.= "</map>";

//		$content.= '<div class="EG_CollapsibleMap_TitleDiv">
//					<div>
//						<a class="EG_CollapsibleMap_TitleAhref" href="javascript:void(0)" onclick="toggleBlockDiv(\'gene_map\', \'gene_CollapsibleMap\');">
//						<img id="gene_CollapsibleMap"  class ="cl_CollapsibleArea_collapsing" src="sites/all/images/empty.gif">
//						<span class="EG_CollapsibleMap_Title">Gene Map</span></a>
//						
//						<div class="EG_CollapsibleMap_HrDiv"> <hr class="EG_CollapsibleMap_Hr"></div>
//					</div>';
//	   $content.= '<div id="gene_map"><img SRC="'.$tmpfname.'" usemap="#mapgenemap" />';
//		$content.='</div>';
		$content .= '<img style="display:block;" SRC="'.file_create_url($tmpfname).'" usemap="#mapgenemap" />';
	return $content;	
	
}
