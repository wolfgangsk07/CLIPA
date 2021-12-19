<?php
//记录IP
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
try{
	if(file_exists("../others/plotips.txt")){
		$iprecord="../others/plotips.txt";
		$handle=fopen($iprecord,"a+"); 
		$str=fwrite($handle,ipaddress()."\r\n");
		fclose($handle);
	}else{
		$iprecord="../others/plotips.txt";
		$handle=fopen($iprecord,"w"); 
		$str=fwrite($handle,time()."\r\n");
		fclose($handle);
	}
}catch(Exception $e){

}

?>