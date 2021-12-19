<?php
	$output=[];
	$files=["../js/functions.js","../js/plot.js","../r/L_functions.R","../r/W_functions.R","../css/index.css"];
	date_default_timezone_set('PRC');
	foreach($files as $v){
		$a=filemtime($v);
		if(time()-$a<30){
			array_push($output,$v."刚刚修改：".date("Y-m-d H:i:s",$a));
		}
	}
	echo(json_encode($output));
?>