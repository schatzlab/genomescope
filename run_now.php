<html>
    <body>
        
<?php

    if( !isset($_POST['code']) ) { echo shell_exec('echo ERROR: No code passed to run.php >> user_data/ERRORS/run.log');}
    $code=$_POST["code"];
    if( !isset($_POST['kmer_length']) ) { echo shell_exec('echo ERROR: No kmer_length passed to run_now.php >> user_data/$code/run.log');}
    if( !isset($_POST['read_length']) ) { echo shell_exec('echo ERROR: No read_length passed to run_now.php >> user_data/$code/run.log');}
    $kmer_length = $_POST["kmer_length"];
    $read_length = $_POST["read_length"];
    $url="analysis.php?code=$code";
    $filename="user_uploads/$code";
    $oldmask = umask(0);
    mkdir("user_data/$code");
    umask($oldmask);
    
    echo shell_exec("Rscript histfitdup_bothplots.R $filename $kmer_length $read_length user_data/$code &> user_data/$code/errors.log &"); 

    header('Location: '.$url);
?>
    </body>
</html>

<!-- <form name="input_code_form" action="run.php" id="analysis_form" method="post"> -->
