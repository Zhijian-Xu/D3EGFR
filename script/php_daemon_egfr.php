<!--------------------------------------------------
            E-mail: xbzhang@simm.ac.cn
            Created: 2015-5-20
            Last Modify: 2015-5-20
            Version 0.01
-------------2015------Xinben-Zhang--------------->

<!-----

---->

<?php
error_reporting(E_ALL);
ini_set('display_errors', TRUE);
ini_set('display_startup_errors', TRUE);
ini_set("memory_limit","-1");

echo "Welcome to D3EGFR Daemon \n";
date_default_timezone_set("Asia/Shanghai");

// Import the configuration file

// Connect to the database
define("RunningFolderPath", "/home/yqyang/D3EGFR/script/");
define("ScriptsFolderPath","/home/yqyang/D3EGFR/script/");

require 'Client.php';
$redis = new Credis_Client('106.14.188.60',6380);
$redis->auth("SariP@ssw0rddddc1201");
$redis->select(1);
			
// Fetch the file information
// 从数据库中读取任务信息

//restart jobs
$job_running_error = $redis->sdiff('D3EGFR_running','D3EGFR_finished','D3EGFR_error');
$job_submited=$redis->lrange('D3EGFR_submited','0','-1');
$job_restart = array_diff($job_running_error,$job_submited);
for($i=0;$i<count($job_restart);$i++){
    if(!strlen($job_restart[$i])>0)
        break;
    echo $job_restart[$i]." restart\n";
    $redis->lpush('D3EGFR_submited',$job_restart[$i]);
    $redis->srem('D3EGFR_running',$job_restart[$i]);
    $redis->del($job_restart[$i]."_finishedtime");
    $redis->del($job_restart[$i]."_begintime");
}

do {
    // 正在等待的任务数
    $jobs_waiting = $redis->llen('D3EGFR_submited');
    while($jobs_waiting>0){

        echo "jobs waiting:" . $jobs_waiting."\n";

		// 提取排在第一位的任务id
        $jobtitle = $redis->rpop('D3EGFR_submited');
        $jobs_waiting = $redis->llen('D3EGFR_submited');
        $redis->sadd('D3EGFR_running',$jobtitle);
        echo "Get job: ".$jobtitle."\n";
        echo "Current Time: " . date(DATE_RFC822) . "\n";
        echo "begin to run job: ".$jobtitle."\n" ;

        $redis->set($jobtitle.'_begintime',date("Y-m-d/H:i:s"));
                
        $mutation=$redis->get($jobtitle.'_mutation');
        $string_command="sh ".ScriptsFolderPath.'/name-modify.sh '.$mutation;
        exec($string_command,$log);
        $mutation=trim(array_shift($log));
        $redis->set($jobtitle.'_mutation',$mutation);
        $result=$redis->exists('D3EGFRAI_'.$mutation);
        if($result>0){
            echo $jobtitle." is finished\n";
            $redis->sadd('D3EGFR_finished',$jobtitle);
            $redis->set($jobtitle.'_finishedtime',date("Y-m-d/H:i:s"));
            continue;
        }

        $string_command="python ".ScriptsFolderPath.'drugResponse_final.py  -m '.$mutation;
        system($string_command);
        if(file_exists(RunningFolderPath."/Prediction/".$mutation.".txt")){
            $file=fopen(RunningFolderPath."/Prediction/".$mutation.".txt","r");
            while(!feof($file)){
                $line=fgets($file);
                $data=explode(" ",trim($line));
                if(count($data)>=7)
                {
                    $redis->set($mutation."_".$data[1]."_score",$data[6]);
                    $redis->set($mutation."_".$data[1]."_response",$data[2]);
                    $redis->set($mutation."_".$data[1]."_crpr",$data[3]);
                    $redis->set($mutation."_".$data[1]."_sd",$data[4]);
                    $redis->set($mutation."_".$data[1]."_pd",$data[5]);
                }
            }
            $redis->sadd("mutationsAi",$mutation);
        }
        else{
            $redis->sadd('D3EGFR_error',$jobtitle);
            $redis->set($jobtitle.'_errormessage',"unrecognized mutation type");
        }
        echo $jobtitle." is finished\n";

        $redis->sadd('D3EGFR_finished',$jobtitle);
        $redis->set($jobtitle.'_finishedtime',date("Y-m-d/H:i:s"));
    }
    try{
        $redis->subscribe("D3EGFR_newjob",function($client,$channel,$message){
                $client->unsubscribe();
        });
    } catch (Exception $e) {
        continue;
    }
}while(1);		

?>
