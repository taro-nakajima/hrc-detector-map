//JavaScript code for simulation of neutron Laue diffraction pattern at HRC
var version = "0.1";
var TOFconst = 2.286;       // TOF at 1 m is 2.286/sqrt(E)

var X0 = 50;
var Y0 = 250;

var length1=200;
var radius=5;

var a_star = new Array(3);
var b_star = new Array(3);
var c_star = new Array(3);

a_star[0]=1.0;
a_star[1]=0.0;
a_star[2]=0.0;


function draw() {
    document.getElementById("verNum").innerHTML=version;
    document.getElementById("verNum2").innerHTML=version;

    draw_DetMap();

}

function draw_DetMap(){


    radius = Number(document.getElementById('a').value);
    a_star[1] = a_star[1]+ Number(document.getElementById('b').value);

    var canvas = document.getElementById('CanvasDetMap');
    var context = canvas.getContext('2d');

    //refresh
    context.clearRect(0, 0, canvas.width, canvas.height);
    context.strokeStyle = "rgb(0, 0, 0)";
    context.lineWidth=1;


    //text 
    context.font = "italic 10px sans-serif";
    context.fillText(a_star[1], X0, Y0+length1);


    // line
    context.strokeStyle = "rgb(255, 0, 0)";
    context.beginPath();
    context.moveTo(X0, Y0);
    context.lineTo(X0+length1, Y0);
    context.stroke();

    // circle
    context.strokeStyle = "rgb(0, 150, 0)";
    var delta=20;
    var limit=5;
    for (let i=0;i<limit;i+=1){
        let PosX=X0+delta*i;
        let PosY=Y0-delta*i**2;
        context.beginPath();
        context.arc(PosX,PosY, radius, 0, 2 * Math.PI);
        context.stroke();
    }


}

