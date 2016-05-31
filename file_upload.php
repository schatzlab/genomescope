<?php
    
    $code = escapeshellcmd($_POST["code_hidden"]);
    $name = "./user_uploads/" . $code;
    
    move_uploaded_file($_FILES['file']['tmp_name'], $name);
    
?>
