<?php
$res=array();
if(file_exists('../others/visitips.txt')){
	$handle = fopen('../others/visitips.txt',"r");//以只读方式打开一个文件
	$i = 0;
	while(!feof($handle)){//函数检测是否已到达文件末尾 
		$tmp=fgets($handle);
		if($i==0){
			$ctime1=date('Y-m-d H:i:s',intval($tmp));
		};
		if($tmp){// 从文件指针中读取一行
			$i++;
		};
	};
	fclose($handle);
	$res['visitips']['times']=$i;
	$res['visitips']['from']=$ctime1;
}else{
	$res['visitips']['times']=0;
	$res['visitips']['from']="Not yet";
};





if(file_exists('../others/plotips.txt')){
	$handle2 = fopen('../others/plotips.txt',"r");//以只读方式打开一个文件
	$j = 0;
	while(!feof($handle2)){//函数检测是否已到达文件末尾 
		
		$tmp=fgets($handle2);
		if($j==0){
			$ctime2=date('Y-m-d H:i:s',intval($tmp));
		};
		if($tmp){// 从文件指针中读取一行
			$j++;
		};
		
	};
	fclose($handle2);
	$res['plotips']['times']=$j;
	$res['plotips']['from']=$ctime2;
	
}else{
	$res['plotips']['times']=0;
	$res['plotips']['from']="Not yet";
}
exit(json_encode($res));
?>