var analysis_path="analysis.php?code=";


//////////////////////////////////////////////////////////////
/////// For analysis page:
//////////////////////////////////////////////////////////////

function showProgress() {
    var run_id_code=getUrlVars()["code"];
    var prog=0;
    
//  remember ajax is asynchronous, so only the stuff inside the success: part will be called after retrieving information. If I put something after the statement, it can't use the info from check_progress.php because it is executed before this php script is called
    //alert('before ajax');
    jQuery.ajax({ 
        type:"POST",
        url: "check_progress.php",
        dataType: 'json',
        data: {code: run_id_code},
        success: function (obj) {
            // alert("inside success");
            // alert(obj);
            prog=obj;

            if (prog=='done') {
                
                //alert("no error in object");
                //
                //alert(prog);
                //alert(prog.length);
                // prog is progress in percent
                
                document.getElementById("plot_info").innerHTML="Analysis done";
                document.getElementById("progress_panel").className = "panel panel-success center";
                check_plot_exists(0);
            }
            else if (prog=="in progress") {
                document.getElementById("plot_info").innerHTML="In progress";
                setTimeout(function() {showProgress()},100);
            }
            else if (prog.indexOf("not converge") > -1) {
                document.getElementById("plot_info").innerHTML=prog + " Try uploading your file and running it again.";
                document.getElementById("progress_panel").className = "panel panel-danger center";
            }
            else {
                //alert("ERROR in getting json data back from check_progress.php to copycat_analysis.js");
                document.getElementById("plot_info").innerHTML = "Running...";
                console.log("else");
                setTimeout(function(){showProgress();},500);
                
            }
        }
        
    });

}



function getUrlVars() {
    var vars = {};
    var parts = window.location.href.replace(/[?&]+([^=&]+)=([^&]*)/gi, function(m,key,value) {
        vars[key] = value;
    });
    return vars;
}
//
//function test() {
//    var run_id_code=getUrlVars()["code"];
//    alert(run_id_code);
//}
//


function check_plot_exists(counter) {
    
    var run_id_code=getUrlVars()["code"];
    var plot_url="user_data/"+run_id_code + "/plot.png";
    
    
    if (counter>100) {
        alert("Taking too long to find "+ plot_url + " counter: " + counter);
    }
    else {
        
        jQuery.ajax({ 
            url: plot_url,
            error: function() {
                console.log(counter+1);
                setTimeout(function(){check_plot_exists(counter+1);},500);

            },
            success: function () {
                document.getElementById("results").style.visibility= 'visible';
                // alert("inside success");
                document.getElementById("landing_for_plot1").innerHTML='<img class="fluidimage" src="' + "./user_data/"+ run_id_code + "/plot.png"  + ' "/>'; // width="750" height="500" // style="max-height:100%; max-width:100%;" 
                document.getElementById("landing_for_plot1_details").innerHTML='<iframe  width="600" height="300" src="' + "./user_data/"+ run_id_code + "/summary.html" + '" frameborder="0"></iframe>';

                document.getElementById("down_txt_1").href = "./user_data/"+ run_id_code + "/summary.txt";
                document.getElementById("down_img_1").href = "./user_data/"+ run_id_code + "/plot.png";
                imageresize();
            }
        });
    }
}


//////////////////////////////////

// $(window).resize(function() {
//     // if(this.resizeTO) clearTimeout(this.resizeTO);
//     // this.resizeTO = setTimeout(function() {
//     //     $(this).trigger('resizeEnd');
//     // }, 500);
//     console.log("resizing")
// });

// //redraw graph when window resize is completed  
// $(window).on('resizeEnd', function() {
//     console.log("resize end")
//     document.getElementById("landing_for_plot1").height = $(window).height;
// });
function imageresize() {
    console.log("resizing")
    var top_padding = 150;
    var side_padding = 0.05;
    var aspect_ratio = 1.5;
    var height = Math.min($( window ).width()/aspect_ratio*(1-side_padding), $( window ).height()-top_padding);
    $(".fluidimage").height(height + "px");
    $(".fluidimage").width(height*aspect_ratio + "px");
}


$(document).ready(function() {
    showProgress();
    

    $(window).bind("resize", function(){//Adjusts image when browser resized
       imageresize();
    });
    
});





// How to execute code after getting info from multiple files:
    //
    //$.when(
    //    $.get(filename_input, function(csvString) {
    //        array_input = $.csv.toArrays(csvString, {onParseValue: $.csv.hooks.castToScalar}); 
    //    }),
    //    $.get(filename_output, function(csvString) {
    //        array_ouput = $.csv.toArrays(csvString, {onParseValue: $.csv.hooks.castToScalar} ); 
    //    })
    //).then(function() {
    //    console.log(array_input)
    //    console.log(typeof array_input)
    //    console.log(array_input.length)
    //    console.log(array_input[0].length)
    //    var diff=[]
    //    for (i=0; i<v_tick_max+2; i++){
    //        v_ticks.push(i);
    //    }
    //});