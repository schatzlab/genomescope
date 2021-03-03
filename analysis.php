<!DOCTYPE html>

<html>


<!--    NAVIGATION BAR-->
    <?php include "header.html";?>
    <link rel="stylesheet" type="text/css" media="screen" href="http://cdnjs.cloudflare.com/ajax/libs/fancybox/1.3.4/jquery.fancybox-1.3.4.css" />

            <!-- ////////////////////////////////////////////////// -->
            <!-- ////////////////      RESULTS     //////////////// -->
            <!-- ////////////////////////////////////////////////// -->
            <div class="thumbnail frame">
                <div class = "caption" style="text-align: center"><h3 id="description_header"></h3></div>
            </div>

            <div id="results">
               
                <!--  All plots  -->
                <div class="thumbnail plot_frame frame">
                    <div style="display:inline-block" id="landing_for_plot1" class="plot_img">
                        <!-- Landing spot for plot image -->
                    </div>
                    <div style="display:inline-block" id="landing_for_plot2" class="plot_img">
                        <!-- Landing spot for plot image -->
                    </div>
                </div>
                <div class="thumbnail plot_frame frame">
                    <div style="display:inline-block" id="landing_for_plot3" class="plot_img">
                        <!-- Landing spot for plot image -->
                    </div>
                    <div style="display:inline-block" id="landing_for_plot4" class="plot_img">
                        <!-- Landing spot for plot image -->
                    </div>
                </div>


                <div class="row">
                    <div class="col-lg-6"> 

                        <!-- Assembly statistics N50 etc. -->
                        <div class="thumbnail frame plot_frame">
                            <div class="caption">
                                <h4>Results</h4>
                                <p id="landing_for_text1"></p>
                            </div>
                        </div>
                    </div>
                    <div class="col-lg-6"> 
                        <!-- Variant statistics -->
                        <div class="thumbnail frame plot_frame">
                            <div class="caption">
                                <h4>Model</h4>
                                <p id="landing_for_text2"></p>
                            </div>
                        </div>
                    </div>
                </div>
                
            </div>
   
    <!--View analysis later-->
    <div id="codepanel" >
        <div class="panel panel-info">
          <div class="panel-heading">
            <h3 class="panel-title">View analysis later</h3>
          </div>
          <div id="code" class="panel-body">
            <?php
                $code=$_GET["code"];
                $url="http://genomescope.org/genomescope2.0/analysis.php?code=$code";
    
                echo "Return to view your results at any time: <input type=\"text\" class=\"form-control\" value=\"$url\"></input>";
            ?>
          </div>
        </div>
    </div>

    <!-- ////////////////////////////////////////////////// -->
    <!-- /////////////      Progress info     ///////////// -->
    <!-- ////////////////////////////////////////////////// -->
    <div class="panel panel-info center" id="progress_panel">
      <div class="panel-heading">
        <h3 class="panel-title">Progress</h3>
      </div>
      <div class="panel-body">
        <div id="plot_info">
        Checking progress...
        </div>
      </div>
    </div> <!-- End of progress info -->
    <!-- ////////////////////////////////////////////////// -->
    


    
<!--   jquery must be first because bootstrap depends on it   -->
<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
<script src="js/bootstrap.min.js"></script>


<script src="js/jquery.csv-0.71.min.js"></script>

<script type="text/javascript" src="https://www.google.com/jsapi"></script>
<script type='text/javascript' src="js/analysis_page_script.js?rndstr="<?php rand(100000,999999) ?> ></script>

<script type="text/javascript" src="http://code.jquery.com/jquery-migrate-1.2.1.min.js"></script>
<script type="text/javascript" src="http://cdnjs.cloudflare.com/ajax/libs/fancybox/1.3.4/jquery.fancybox-1.3.4.pack.min.js"></script>

</body>
</html>




