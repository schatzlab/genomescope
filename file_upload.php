<?php
    
    
    $code = $_POST["code_hidden"];
    $name = "./user_uploads/" . $code;
    
    move_uploaded_file($_FILES['file']['tmp_name'], $name);

    ////for debugging:
    //file_put_contents( 'yowtf', print_r($_POST["code_hidden"], true));
    //
    //
    //file_put_contents( 'yohai', print_r($name, true));
    //
    
?>