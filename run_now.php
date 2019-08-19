<html>
    <body>
        
<?php

    if( !isset($_POST['code']) ) { echo shell_exec('echo ERROR: No code passed to run.php >> user_data/ERRORS/input_validation.log');}
    $code=$_POST["code"];
    mkdir("user_data/$code");

    if( !isset($_POST['kmer_length']) ) { echo shell_exec("echo ERROR: No kmer_length passed to run_now.php >> user_data/$code/input_validation.log");}
    if( !isset($_POST['ploidy']) ) { echo shell_exec("echo ERROR: No ploidy passed to run_now.php >> user_data/$code/input_validation.log");}
    if( !isset($_POST['max_kmer_cov']) ) { echo shell_exec("echo ERROR: No max_kmer_cov passed to run_now.php >> user_data/$code/input_validation.log");}
    if( !isset($_POST['lambda']) ) { echo shell_exec("echo ERROR: No lambda passed to run_now.php >> user_data/$code/input_validation.log");}
    if( !isset($_POST['description']) ) { echo shell_exec("echo ERROR: No description passed to run_now.php >> user_data/$code/input_validation.log");}

    echo shell_exec("echo user_data/$code/input_validation.log >> user_data/$code/input_validation.log");

    $kmer_length = escapeshellarg($_POST["kmer_length"]);
    $ploidy = escapeshellarg($_POST["ploidy"]);
    $max_kmer_cov = escapeshellarg($_POST["max_kmer_cov"]);
    $lambda = escapeshellarg($_POST["lambda"]);
    $description = escapeshellarg($_POST["description"]);

    echo shell_exec("echo $description > user_data/$code/description.txt");

    $url="analysis.php?code=$code";
    $filename="user_uploads/$code";
    $oldmask = umask(0);
    umask($oldmask);

    // For old version:
    // echo shell_exec("Rscript histfitdup_bothplots.R $filename $kmer_length $read_length user_data/$code &> user_data/$code/errors.log &"); 
    
    // For new version:
    // genomescope.R histogram_file k-mer_length read_length output_dir
    // echo shell_exec("Rscript genomescope.R $filename $kmer_length $read_length user_data/$code $max_kmer_cov &> user_data/$code/run.log &"); 

    // For 2.0:
    // genomescope.R -i histogram_file -k k-mer_length -p ploidy -o output_dir
    echo shell_exec("./genomescope.R -i $filename -k $kmer_length -p $ploidy -l $lambda -o user_data/$code &> user_data/$code/run.log &");

    $new_dataset = array( "date"=>time(), "codename"=>$code, "description"=> $description );

    $my_datasets = array();
    if(isset($_COOKIE["results"])) {
      // echo "cookie is already there, adding to it.";
      $my_datasets = json_decode($_COOKIE["results"], true);
    } else {
      // echo "cookie not set, creating new one";
    }
    array_push($my_datasets, $new_dataset);
    setcookie("results", json_encode($my_datasets));


    header('Location: '.$url);
?>
    </body>
</html>

<!-- <form name="input_code_form" action="run.php" id="analysis_form" method="post"> -->
