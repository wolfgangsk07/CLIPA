<?php
function ipaddress() {
                global $ip;
                if(getenv("HTTP_CLIENT_IP"))
                        $ip=getenv("HTTP_CLIENT_IP");
                else if(getenv("HTTP_X_FORWARDED_FOR"))
                        $ip=getenv("HTTP_X_FORWARDED_FOR");
                else if(getenv("REMOTE_ADDR"))
                        $ip=getenv("REMOTE_ADDR");
                else
                        $ip="Unknow";
               return $ip;
}
try
{
	list($msec, $sec) = explode(' ', microtime());
	$msectime = (float)sprintf('%.0f', (floatval($msec) + floatval($sec)) * 1000);

	//set_time_limit(10);
	$userParams=base64_decode($_POST['userParams']);
	$userParams_obj=json_decode($userParams);

	$use=$userParams_obj -> use;
	$requestGUID=$userParams_obj -> requestGUID;
	$path=$userParams_obj -> path;
	$userParams_path="../plotCache/userParams_".$requestGUID.".json";
	$up = fopen($userParams_path,"w");
		fwrite($up,$userParams);
		fclose($up);
	$tasklist=array();
	if(file_exists("c:/windows")){
		exec("tasklist",$tasklist);
	}else{
		exec("ps -aux",$tasklist);
	}
	$runningcount=0;
	for($t=0;$t<count($tasklist);$t++){
		if(strpos($tasklist[$t],"Rscript")!==false){
			$runningcount++;
		}
	}
	if($runningcount>8){
		$res=array();
		$res['requestGUID']=$requestGUID;
		$res['res']['error']=$runningcount." task(s) are running currently, please try again later.";
		exit(json_encode($res));
	};
	
	$shell="Rscript ../r/centralHub.r ".$requestGUID;
	
	$pid=exec($shell);
	//echo $pid;

	for($i=0;$i<3;$i++){
		if(file_exists($path)){
			sleep(0.1);
			$fp = fopen($path,"r");
				$str = fread($fp,filesize($path));//指定读取大小，这里把整个文件内容读取出来
				echo $str;
				fclose($fp);
				break;
		}
		sleep(1);
	}
	
}catch(Exception $e){
    echo 'Message: ' .$e;
}

?>