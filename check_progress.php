<?php
    $code=$_POST['code'];
    
    $filename="user_data/" . $code . "/progress.txt";
    
    if (file_exists($filename)) {
        $myfile = fopen($filename, "r") or die("error");
        
        $progress_stats = "in progress";
        
        while(!feof($myfile)) {
            $line=fgets($myfile);
            $line=trim(preg_replace( '/\s+/', ' ', $line ));
            if (strlen($line)>0) {
                $progress_stats = $line;
            }
        }
        fclose($myfile);
        echo json_encode($progress_stats);
    }
    else {
        $progress_stats = "file " + $filename + " doesn't exist";
        echo json_encode($progress_stats);
    }

?>
