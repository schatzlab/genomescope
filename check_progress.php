<?php
    $code=$_POST['code'];
    
    $filename="user_data/" . $code . "/progress.txt";
    
    if (file_exists($filename)) {
        $progress_stats = file($filename);
        echo json_encode($progress_stats);
    }
    else {
        $progress_stats = array("File " . $filename . " doesn't exist.");
        echo json_encode($progress_stats);
    }

?>
