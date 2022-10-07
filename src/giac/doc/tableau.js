/**
 * This class represents a polyline drawn by a user
 */
class Delineation {
  color="";
  points = [];
  lastpoints = [];
  pressures=[];
  lastpressure='';
  
  reset() {
    this.drawing = true;
    this.lastpoints = this.points;
    this.points = [];
    this.pressures = [];
    this.lastpressure='';
  }

  finish() {
    //console.trace();
    this.drawing = false;
    this.removePolygon();
  }
  /**
   * @returns true if the set of current points is non-empty
   */
  isDrawing() {
    return this.points.length > 0;
  }

  containsPolygonToMagnetize() {
    return this.points.length > 0;
  }

  drawPolygon(points) { 
    if (document.getElementById("magnetCreationPolygon"))
      return;

    let polyline = document.createElementNS('http://www.w3.org/2000/svg', 'polyline');
    polyline.id = "magnetCreationPolygon";
    svg.appendChild(polyline);
    console.log('drawPolygon',svg);

    points.push(points[0]);
    polyline.setAttribute("points", points.map((p) => p.x + ", " + p.y).join(" "));
  }

  removePolygon() {
    if (document.getElementById("magnetCreationPolygon")){
      console.log('removePolygon',svg); //console.trace();
      svg.removeChild(document.getElementById("magnetCreationPolygon"));
    }
  }


  addPoint(point,pressure=0.5) {
    let l=this.points.length,lastp=l?this.points[l-1]:0;
    // let dotinprevpoly=this.isDot() && this.dotInPreviousPolygon();
    let dotinprevpoly= this.isDot() && Delineation.inPolygon(point,this.lastpoints);
    if (l && !dotinprevpoly && point.x==lastp.x && point.y==lastp.y)
      return;
    if (pressure=='')
      pressure=0.3;
    let curp=Math.floor(pressure*10-0.0001);
    if (curp<0) curp=0;
    if (curp>9) curp=9;
    // console.log(point.x,point.y,lastp.x,lastp.y);
    if (l>=2 && !dotinprevpoly){
      let p1=this.points[l-2];
      //console.log(p1.x,p1.y,lastp.x,lastp.y,point.x,point.y);
      if (UI.p2inp1p3(p1.x,p1.y,lastp.x,lastp.y,point.x,point.y,0.3) && curp==this.lastpressure){
	this.points[l-1]=point;
	return;
      }
    }
    this.points.push(point);
    this.pressures.push(curp);
    this.lastpressure=curp;
    if (dotinprevpoly) {
      this.drawPolygon(this.lastpoints);

      window.setTimeout(() => {
	// console.log('dotinprevpoly 2',this.drawing,this.isDot(),this.dotInPreviousPolygon());
        if (this.drawing &&
	    this.isDot() && this.dotInPreviousPolygon()) {
          this.removePolygon();
          this._removeContour(); //remove the dot
          this.points = this.lastpoints;
          this.lastpoints = [];
	  //BoardManager.cancelStack.stack.pop();
          this.cutAndMagnetize();
        }
      }, 2000);
    }
    else
      this.removePolygon();

  }

  /**
   * @returns true if the current drawing is just a point
   */
  isDot() {
    if (this.points.length == 0)
      return false;

    for (let point of this.points)
      if (Math.abs(point.x - this.points[0].x) > 1 && Math.abs(point.y - this.points[0].y) > 1)
        return false;

    return true;
  }

  /**
   * 
   * @param {*} point 
   * @param {*} polygon 
   * @returns true iff the point is inside the polygon
   */
  static inPolygon(point, polygon) {
    // ray-casting algorithm based on
    // https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html/pnpoly.html
    if (!point) return false;
    var x = point.x, y = point.y;

    var inside = false;
    for (var i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
      var xi = polygon[i].x, yi = polygon[i].y;
      var xj = polygon[j].x, yj = polygon[j].y;

      var intersect = ((yi > y) != (yj > y))
          && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
      if (intersect) inside = !inside;
    }

    return inside;
  }

  dotInPreviousPolygon() {
    return Delineation.inPolygon(this.points[0], this.lastpoints);
  }

  erase() {
    if (!this.isSuitable())
      return;

    this._removeContour();
    this._clearBehindMagnet();
    BoardManager.save();
  }


  /**
   * @description magnetize the "selected" part of the blackboard. The selected part is also removed.
   */
  cutAndMagnetize() {
    if (!this.isSuitable())
      return;
    console.log('cutmagnetize',this);
    BoardManager.cancelStack.stack.pop();
    this._removeContour();
    this._createMagnetFromImg();
    this._clearBehindMagnet();
    //BoardManager.cancelStack.stack.pop();
    BoardManager.save();
  }


  /**
   * @description magnetize the "selected" part of the blackboard.
   */
  copyAndMagnetize() {
    if (!this.isSuitable())
      return;
    console.log('copymagnetize',this);

    this._removeContour();
    this._createMagnetFromImg();
    BoardManager.save();
  }




  isSuitable() {
    for (let point of this.points) {
      if (Math.abs(point.x - this.points[0].x) > 16 &&
          Math.abs(point.x - this.points[0].x) > 16)
        return true;
    }
    return false;
  }

  _removeContour = () => {
    const context = document.getElementById("canvas").getContext("2d");
    context.globalCompositeOperation = "destination-out";
    context.strokeStyle = "rgba(255, 255, 255, 1)";
    context.lineWidth = 6;
    context.globalAlpha = 1.0;


    context.moveTo(this.points[0].x, this.points[0].y);
    for (let point of this.points) {
      context.lineTo(point.x, point.y);
    }
    context.stroke();
  }

  _getRectangle() {
    let r = { x1: canvas.width, y1: canvas.height, x2: 0, y2: 0 };

    for (let point of this.points) {
      r.x1 = Math.min(r.x1, point.x);
      r.y1 = Math.min(r.y1, point.y);
      r.x2 = Math.max(r.x2, point.x);
      r.y2 = Math.max(r.y2, point.y);
    }

    return r;
  }




  _createMagnetFromImg = () => {
    const rectangle = this._getRectangle();
    console.log('createmagnetfromimg',rectangle);
    if (1){ 
      //UI.console_show_stack();
      let xmin=rectangle.x1,xmax=rectangle.x2,ymin=rectangle.y1,ymax=rectangle.y2;
      let delineations=UI.cut_stack(xmin,xmax,ymin,ymax,this.points);
      console.log(delineations);
      if (delineations!=''){
	// BoardManager.cancelStack.stack=tab; // needed?
	//BoardManager.cancelStack.back();
	//UI.console_show_stack();
	BoardManager.cancelStack.back();
	BoardManager.redraw_handwritings();
	//BoardManager.cancelStack.currentIndex=tab.length;
	//console.log(delineations);
	MagnetManager.addMagnetText(xmin, ymin,'','handwriting',delineations,0);
	return;
      }
    }
    let img = new Image();
    BoardManager._toBlobOfRectangle(rectangle, (blob) => img.src = URL.createObjectURL(blob));
    //img.src = BoardManager._getDataURLPictureOfRectangle(rectangle);
    img.style.clipPath = "polygon(" + this.points.map(point => `${point.x - rectangle.x1}px ${point.y - rectangle.y1}px`).join(", ") + ")";
    MagnetManager.addMagnet(img);
    img.style.left = rectangle.x1 + "px";
    img.style.top = rectangle.y1 + "px";
  }

  _clearBehindMagnet = () => {
    const context = document.getElementById("canvas").getContext("2d");
    context.save();
    context.beginPath();
    context.moveTo(this.points[0].x, this.points[0].y);
    for (let point of this.points) {
      context.lineTo(point.x, point.y);
    }
    context.clip();
    context.clearRect(0, 0, window.innerWidth, window.innerHeight);
    context.restore();
    context.globalCompositeOperation = "source-over";
    this.reset();
  }



}

// Palette.js
class Palette {
  /** colors that can have a chalk. The first color must be white */
  colors = ["white", "yellow", "orange", "rgb(100, 172, 255)", "Crimson", "Plum", "LimeGreen"];

  buttons = [];
  currentColorID = 0;
  onchange = () => { };
  

  static radius = 96;

  _createPalette() {
    const div = document.getElementById("palette");
    for (let i in this.colors) {
      this.buttons[i] = this._createColorButton(i);
      div.appendChild(this.buttons[i]);
    }
  }

  /**
   * @description switch the first color (white <=> black)
   */
  switchBlackAndWhite() {
    this.colors[0] = (this.colors[0] == "white") ? "black" : "white";
    this.onchange();
    UI.caseval_noautosimp('color('+this.colors[0]+')');
  }

  /**
   * 
   * @param {*} i  an index between 0 and this.colors.length - 1
   * @description create the button for the color of index i
   */
  _createColorButton(i) {
    const img = new Image();
    img.src = ChalkCursor.getCursorURL(this.colors[i]);
    img.classList.add("paletteColorButton");

    let angle = -Math.PI / 2 + 2 * Math.PI * i / this.colors.length;

    img.style.top = (Palette.radius * Math.sin(angle) - 22) + "px";
    img.style.left = (Palette.radius * Math.cos(angle) - 16) + "px";
    img.style.borderColor = this.colors[i];

    img.onclick = () => {
      this.buttons[this.currentColorID].classList.remove("selected");
      this.currentColorID = i;
      this.buttons[this.currentColorID].classList.add("selected");
      this.hide();
      this.onchange();
    }
    return img;
  }

  /**
   * @description select the next color
   */
  next() {
    this.buttons[this.currentColorID].classList.remove("selected");
    this.currentColorID++;
    this.currentColorID = this.currentColorID % this.colors.length;
    this.buttons[this.currentColorID].classList.add("selected");
    this.onchange();
  }

  /**
   * @description select the previous color
   */
  previous() {
    this.buttons[this.currentColorID].classList.remove("selected");
    this.currentColorID--;
    if (this.currentColorID < 0)
      this.currentColorID = this.colors.length - 1;
    this.buttons[this.currentColorID].classList.add("selected");
    this.onchange();
  }

  /**
   * @param position a point {x: ..., y: ...}
   * @description show the palette at position position
   */
  show(position) {
    const div = document.getElementById("palette");
    div.innerHTML = "";
    this._createPalette();

    position.y = Math.max(position.y, Palette.radius + 16 + 48);
    position.x = Math.max(position.x, Palette.radius + 16 + 48);

    div.style.left = position.x + "px";
    div.style.top = position.y + "px";
    div.classList.remove("PaletteHide");
    div.classList.add("PaletteShow");
  }

  hide() {
    const div = document.getElementById("palette");
    div.classList.remove("PaletteShow");
    div.classList.add("PaletteHide");
  }

  isShown() {
    return document.getElementById("palette").classList.contains("PaletteShow");
  }

  /**
   * @returns the selected color
   */
  getCurrentColor() {
    return this.colors[this.currentColorID];
  }
}

const ERASEMODEDEFAULTSIZE = 10;

class User {
  xInit = 0;
  yInit = 0;

  x = 0;
  y = 0;
  isDrawing = false;
  alreadyDrawnSth = false; // true if something visible has been drawn (If still false, draw a dot)
  eraseMode = false;
  eraseModeBig = false;
  lastDelineation = new Delineation();
  canWrite = true;

  color = "white";

  cursor = undefined;

  userID = 0;

  setUserID(userID) {
    this.userID = userID;
  }

  setCanWrite(bool) {
    this.canWrite = bool;
  }

  /**
   * 
   * @param {*} isCurrentUser that tells whether the user is the current one
   * @description create the user. 
   */
  constructor(isCurrentUser) {
    this.cursor = document.createElement("div");
    this.cursor.classList.add("cursor");

    if (isCurrentUser)
      this.cursor.hidden = true;

    document.getElementById("cursors").appendChild(this.cursor);
    if (isCurrentUser)
      document.getElementById("canvas").style.cursor = ChalkCursor.getStyleCursor(this.color);
  }



  updateCursor() {
    if (this.isCurrentUser()) {
      document.getElementById("canvas").style.cursor = ChalkCursor.getStyleCursor(this.color);
    }
  }

  /**
   * tells that the user has disconnected
   */
  destroy() {
    document.getElementById("cursors").removeChild(this.cursor);
  }

  setCurrentColor(color) {
    //console.log('setcurrentcolor',color);
    this.color = color;
    this.updateCursor();
  }

  getCurrentColor() {
    return this.color;
  }

  switchChalkEraser() {
    this.eraseMode = !this.eraseMode;

    if (this.eraseMode) {

    }
    else {

    }
  }


  switchChalk() {
    this.eraseMode = false;

    if (this.isCurrentUser()) {
      this.updateCursor();
      buttonEraser.innerHTML = "⌫";
    }

  }

  /**
   * @returns true iff the user is the current user (the one that controls the mouse)
   */
  isCurrentUser() {
    return (this == user);
  }

  switchErase() {
    this.eraseMode = true;

    if (this.isCurrentUser()) {
      palette.hide();
      document.getElementById("canvas").style.cursor = EraserCursor.getStyleCursor();
      buttonEraser.innerHTML = UI.langue==-1?"Craie":"Chalk";
    }
  }


  mousedown(evt) {
    MagnetManager.setInteractable(false);

    //unselect the selected element (e.g. a text in edit mode)
    document.activeElement.blur();


    //console.log("mousedown")
    this.x = evt.offsetX;
    this.y = evt.offsetY;
    this.xInit = this.x;
    this.yInit = this.y;
    this.isDrawing = true;
    this.eraseModeBig = false;

    if (this.canWrite) {
      if (this.eraseMode) {
	let delineations=UI.cut_stack(this.x-ERASEMODEDEFAULTSIZE/2,this.x+ERASEMODEDEFAULTSIZE/2,this.y-ERASEMODEDEFAULTSIZE/2,this.y+ERASEMODEDEFAULTSIZE/2);
	if (delineations!='')
	// console.log('cut',delineations,this.x,this.y,ERASEMODEDEFAULTSIZE);
	  BoardManager.redraw_handwritings();
        // clearLine(this.x, this.y, this.x, this.y, ERASEMODEDEFAULTSIZE);
      }
      else {
        this.lastDelineation.reset();
	this.lastDelineation.color=user.getCurrentColor(); 
        this.lastDelineation.addPoint({ x: this.x, y: this.y });
      }

    }


    palette.hide();
  }



  mousemove(evt) {
    let evtX = evt.offsetX;
    let evtY = evt.offsetY;
    let pressure=evt.pressure;
    if (pressure=='' || pressure==0)
      pressure=0.3;
    this.cursor.style.left = evtX - 8;
    this.cursor.style.top = evtY - 8;

    if (this.canWrite) {
      if (this.isDrawing) {//} && this.lastDelineation.isDrawing()) {
        palette.hide();
        if (this.eraseMode) {
          let lineWidth = 10;

          lineWidth = 10 + 30 * pressure;

          if (Math.abs(this.x - this.xInit) > window.innerWidth / 4 || Math.abs(this.y - this.yInit) > window.innerHeight / 4)
            this.eraseModeBig = true;

          if (this.eraseModeBig) {
            lineWidth = 128;
          }

          if (this.isCurrentUser())
            document.getElementById("canvas").style.cursor = EraserCursor.getStyleCursor(lineWidth);
          let delineations=UI.cut_stack(this.x-lineWidth/2,this.x+lineWidth/2,this.y-lineWidth/2,this.y+lineWidth/2);
	  if (delineations!='') BoardManager.redraw_handwritings();
        }
        else {

          if (this.lastDelineation.isDrawing()) {//this guard is because, when a magnet is created the user does not know the drawing stopped.
            drawLine(document.getElementById("canvas").getContext("2d"), this.x, this.y, evtX, evtY, pressure, this.color);
            this.lastDelineation.addPoint({ x: evtX, y: evtY },pressure);
	    //console.log('move',evtX,evtY,pressure);
          } 
          

          

        }

        if (Math.abs(this.x - this.xInit) > 1 || Math.abs(this.y - this.yInit) > 1)
          this.alreadyDrawnSth = true;
      }
    }

    this.x = evtX;
    this.y = evtY;
  }


  mouseup(evt) {
    MagnetManager.setInteractable(true);
    // alert(evt.pressure);
    if (this.canWrite) {
      this.lastDelineation.finish();

      if (this.lastDelineation.points.length){
	let colors=this.lastDelineation.color;
	if (colors=="white")
	  colors='';
	else
	  colors +=",";
	//console.log('mouseup',this.lastDelineation.color,colors);
	BoardManager.cancelStack.push(colors+this.lastDelineation.points[0].x+","+this.lastDelineation.points[0].y+","+UI.encode(this.lastDelineation.points,this.lastDelineation.pressures));
      }
      UI.link(0);
      //console.log("mouseup",this.lastDelineation.points);
      //console.log("mouseup",this.lastDelineation.color+","+this.lastDelineation.points[0].x+","+this.lastDelineation.points[0].y+","+UI.encode(this.lastDelineation.points,this.lastDelineation.pressures));
      if (this.isDrawing && !this.eraseMode && !this.alreadyDrawnSth) {
        drawDot(this.x, this.y, this.color);
      }

      if (this.isCurrentUser()) {
        if (this.eraseMode) //restore the eraser to the original size
          document.getElementById("canvas").style.cursor = EraserCursor.getStyleCursor();
      }



      BoardManager.saveCurrentScreen();
    }
    this.alreadyDrawnSth = false;
    this.isDrawing = false;
  }
}

window.onload = load;
window.onresize = resize;

let users = {};
let user = undefined;
let palette = new Palette();

/**
 * this function sets the document.body scrolls to 0
 * It solves some issues:
 * - on smartphones: that we can scroll the page itself
 * - when typing texts in magnet, it makes the screen not to scroll
 */
function installBodyNoScroll() {
	setInterval(() => {
		document.body.scrollLeft = 0;
		document.body.scrollTop = 0;
	}, 1000);
}

function load() {
  installBodyNoScroll();

  user = new User(true);
  users['root'] = user;
  user.setUserID('root');

  buttonNoBackground.onclick = () => { backgroundClear(); Menu.hide(); };
  buttonMusicScore.onclick = () => { musicScore(); Menu.hide(); };

  LoadSave.init();
  BoardManager.init();
  Menu.init();
  Share.init();

  if (navigator.userAgent.match(/Android/i)
      || navigator.userAgent.match(/webOS/i)
      || navigator.userAgent.match(/iPhone/i)
      || navigator.userAgent.match(/iPad/i)
      || navigator.userAgent.match(/iPod/i)
      || navigator.userAgent.match(/BlackBerry/i)
      || navigator.userAgent.match(/Windows Phone/i)) {
    buttonCloseControls.hidden = true; //on phone or tablet, we can not remove the toolbar
  }

  let changeColor = () => {
    if (MagnetManager.getMagnetUnderCursor() == undefined) { //if no magnet under the cursor, change the color of the chalk
      //if (!isDrawing)
      palette.show({ x: user.x, y: user.y });
      palette.next();
    }
    else { // if there is a magnet change the background of the magnet
      let magnet = MagnetManager.getMagnetUnderCursor();
      // console.log(magnet.style.backgroundColor);
      UI.change_magnet_color(magnet);
      UI.link();
    }
  }

  let previousColor = () => {
    if (MagnetManager.getMagnetUnderCursor() == undefined) { //if no magnet under the cursor, change the color of the chalk
      eraseMode = false;

      if (!Delineation.isDrawing)
	palette.show({ x: user.x, y: user.y });
      palette.previous();
    }
    else { // if there is a magnet change the background of the magnet
      let magnet = MagnetManager.getMagnetUnderCursor();
      magnet.style.backgroundColor = previousBackgroundColor(magnet.style.backgroundColor);
    }
  }

  let switchChalkEraser = () => {
    if (user.eraseMode)
      Share.execute("switchChalk", [user.userID]);
    else
      Share.execute("switchErase", [user.userID]);


  }


  //buttonMenu.onclick = Menu.toggle;
  //buttonLeft.onclick = BoardManager.left;
  //buttonRight.onclick = BoardManager.right;
  // buttonCancel.onclick = BoardManager.cancel;
  //buttonMenu2.onclick = Menu.toggle;
  //buttonLeft2.onclick = BoardManager.left;
  //buttonRight2.onclick = BoardManager.right;
  // buttonCancel2.onclick = BoardManager.cancel;
  buttonColors.onclick = changeColor;
  buttonEraser.onclick = switchChalkEraser;
  buttonDivide.onclick = divideScreen;
  //buttonRedo.onclick = BoardManager.redo;

  let buttons = document.getElementById("controls").children;

  for (let i = 0; i < buttons.length; i++)
    buttons[i].onfocus = document.activeElement.blur;

  Welcome.init();


  BlackVSWhiteBoard.init();

  palette.onchange = () => {
    Share.execute("switchChalk", [user.userID]);
    Share.execute("setCurrentColor", [user.userID, palette.getCurrentColor()]);
  }


  document.onkeydown = (evt) => {
    //console.log("ctrl: " + evt.ctrlKey + " shift:" + evt.shiftKey + " key: " + evt.key)
    //if (evt.key == "Backspace") evt.preventDefault();
    if (UI.focused && evt.ctrlKey && evt.key == "m" ){
      $id('assistant_cas').style.display='block';
      evt.preventDefault();
      return;
    }
    if (UI.tableaucm==-1) { UI.tableaucm=0; return; }
    // console.log('ctrl',UI.tableaucm,evt.ctrlKey,evt.shiftKey);
    if (evt.key == "Escape" || evt.key == "F1") {//escape => show menu
      if (Welcome.isShown()){
	Welcome.hide();
	Menu.hide();
      }
      else if (palette.isShown())
	palette.hide();
      else
	Menu.toggle();
      return;
    }
    if (UI.tableaucm || $id('online_help').style.display!='none' || Menu.isShown()) return ;
    if (evt.ctrlKey && evt.key=="s"){
      LoadSave.save();
      return;
    }
    if (evt.key == "h"){  //h=hide toolbar
      welcome.hidden = true;
      controls.hidden=!controls.hidden;
      Welcome.hide();
      Menu.hide();
    }
    if (evt.key == "d" || evt.key=="2"){  //d = divide screen/2
      divideScreen();
      Welcome.hide();
      Menu.hide();
    }
    if (evt.key=="F9"){
      UI.exec_magnets();
      Welcome.hide();
      Menu.hide();
      return;
    }
    if (evt.key=="F8"){
      UI.reorder_magnets();
      Welcome.hide();
      Menu.hide();
      return;
    }
    if (evt.key == "3"){  //3 = divide screen /3
      divideScreen(3);
      Welcome.hide();
      Menu.hide();
    }
    if (evt.key == "Enter" && palette.isShown())
      palette.hide();
    if (evt.key == "ArrowLeft"){
      if (palette.isShown())
	palette.previous();
      else
	BoardManager.left();
    }
    if (evt.key == "ArrowRight"){
      if (palette.isShown())
	palette.next();
      else
	BoardManager.right();
    }
    if (!evt.ctrlKey){
      if (evt.key=="r"){
	UI.set_readonly(true);
	return;
      }
      if (evt.key=="w"){
	UI.set_readonly(false);
	return;
      }
      var evaluator='';
      if (evt.key == "k" || evt.key=="x")  //x = eval with Xcas
	evaluator='xcas';
      if (evt.key == "p")  //p= eval with MicroPython
	evaluator='python';
      if (evt.key == "j")  //j = eval with Javascript
	evaluator='javascript';
      if (evt.key=="x" || evt.key == "k" || evt.key == "p" || evt.key == "j") {
	Welcome.hide();
	UI.ckswitch(evaluator,0);
	// console.log(user.x,user.y); // user.x/user.y are not updated in readonly mode
	UI.addmagnet(user.x,user.y,evaluator);
	evt.preventDefault(); //so that it will not add "new line" in the text element
	return;
      }
    }
    if (evt.key == "l" || evt.key == "$") {
      UI.ckswitch('comment',0);
      Welcome.hide();
      Menu.hide();
      //console.log(user.x,user.y);
      // MagnetManager.addMagnetText(user.x, user.y,'div','latex');
      MagnetManager.addMagnetText(user.x, user.y,'///texte $x^2$','comment',evt.key=="$"?'$':'');
      evt.preventDefault(); //so that it will not add "new line" in the text element
      return;
    }
    if (Menu.isShown()
	// || Welcome.isShown()
       ) return;

    if (!evt.ctrlKey && evt.key == "s"){
      UI.console_show_stack();
    }
    if (!evt.ctrlKey && evt.shiftKey && evt.key == "C") // c => change color
      changeColor();
    if (!evt.ctrlKey && !evt.shiftKey && evt.key == "c"){
      let col=user.getCurrentColor();
      //console.log(col);
      if (col=='white') user.setCurrentColor("yellow");
      if (col=='yellow') user.setCurrentColor("orange");
      if (col=='orange') user.setCurrentColor("white");
      if (col=="Crimson") user.setCurrentColor("Plum");
      if (col=='Plum') user.setCurrentColor("LimeGreen");
      if (col=='LimeGreen') user.setCurrentColor("Crimson");
      // previousColor();
    }
    if (UI.readonly){
      if (evt.key=="y"){
	BoardManager.redo();
	evt.preventDefault();      
      }
      if (evt.key==" "){
	BoardManager.redo(5);
	evt.preventDefault();      
      }
      if (evt.key=="z"){
	BoardManager.cancel();
	evt.preventDefault();      
      }
    }
    if ((evt.ctrlKey && evt.shiftKey && evt.key == "Z") || (evt.ctrlKey && evt.key == "y")) { //ctrl + shift + z OR Ctrl + Y = redo
      BoardManager.redo();
      evt.preventDefault();
    }
    if (evt.ctrlKey && evt.key == "z") {// ctrl + z = undo 
      BoardManager.cancel();
      evt.preventDefault();
      User.lastDelineation.finish();
      // User.lastDelineation.reset();
      
    }
    if (evt.key == "e")  //e = switch eraser and chalk
      switchChalkEraser();
    if (evt.ctrlKey && evt.key == 'x') {//Ctrl + x 
      palette.hide();
      if (user.lastDelineation.containsPolygonToMagnetize()){
	user.lastDelineation.cutAndMagnetize();
	return;
      }
    }
    if (evt.ctrlKey && evt.key == 'c') {//Ctrl + c
      palette.hide();
      if (user.lastDelineation.containsPolygonToMagnetize())
	user.lastDelineation.copyAndMagnetize();
    }
    if (evt.ctrlKey && evt.key == "v") { //Ctrl + v = print the current magnet
      palette.hide();
      MagnetManager.printCurrentMagnet();
    }
    if (evt.key == "m") { //m = make new magnets
      palette.hide();
      if (user.lastDelineation.containsPolygonToMagnetize())
	user.lastDelineation.cutAndMagnetize();
      else {
	MagnetManager.printCurrentMagnet();
	MagnetManager.removeCurrentMagnet();
      }
    }
    if (evt.ctrlKey && evt.key == "p") { //ctrl-p = print the current magnet
      palette.hide();
      MagnetManager.printCurrentMagnet();
    }
    if (evt.key == "Delete" || (evt.ctrlKey && evt.key == "x") || evt.key == "Backspace") { //supr = delete the current magnet
      palette.hide();
      /*if (lastDelineation.containsPolygonToMagnetize())
	lastDelineation.erase();
	else*/
      if (MagnetManager.magnetUnderCursor)
	MagnetManager.currentMagnet=MagnetManager.magnetUnderCursor;
      MagnetManager.removeCurrentMagnet();
      UI.link();
      evt.preventDefault();
    }

  };


  document.getElementById("canvas").onpointerdown = (evt) => {
    if (Welcome.isShown())
      Welcome.hide();
    //console.log('canvas pointerdown',UI.focused);
    if (UI.readonly || UI.focused!=0 ) return;
    evt.preventDefault();
    Share.execute("mousedown", [user.userID, evt])
  };
  document.getElementById("canvasBackground").onpointermove = (evt) => { console.log("mousemove on the background should not occur") };
  document.getElementById("canvas").onpointermove = (evt) => { if (UI.readonly || UI.focused!=0) return; evt.preventDefault(); Share.execute("mousemove", [user.userID, evt]) };
  document.getElementById("canvas").onpointerup = (evt) => { if (UI.readonly || UI.focused!=0) return; evt.preventDefault(); Share.execute("mouseup", [user.userID, evt]) };
  //document.getElementById("canvas").onmousedown = mousedown;

  TouchScreen.addTouchEvents(document.getElementById("canvas"));




  //	document.getElementById("canvas").onmouseleave = function (evt) { isDrawing = false; }

  let magnets = document.getElementsByClassName("magnet");
  for (let i = 0; i < magnets.length; i++) {
    magnets[i].style.top = 0;
    magnets[i].style.left = i * 64 + 2;

  }

  MagnetManager.init();
  loadMagnets();

  BoardManager.load();
  resize();


}


function resize() {
  //alert('resize('+window.innerWidth+','+window.innerHeight+')');
  BoardManager.resize(window.innerWidth, window.innerHeight);
}


function drawLine(context, x1, y1, x2, y2, pressure = 1.0, color = user.getCurrentColor()) {
  context.beginPath();
  context.strokeStyle = color;
  context.globalCompositeOperation = "source-over";
  context.globalAlpha = 0.75 + 0.25 * pressure;
  context.lineWidth = 1 + 3.5 * pressure;
  //console.log(x1,y1,x2,y2,color,context.lineWidth);
  context.moveTo(x1, y1);
  context.lineTo(x2, y2);
  context.stroke();
  context.closePath();
}


function drawDot(x, y, color) {
  const context = document.getElementById("canvas").getContext("2d");
  context.beginPath();
  context.fillStyle = color;
  context.lineWidth = 2.5;
  context.arc(x, y, 2, 0, 2 * Math.PI);
  context.fill();
  context.closePath();
}


function clearLine(x1, y1, x2, y2, lineWidth = 10) {
  const context = document.getElementById("canvas").getContext("2d");
  context.beginPath();
  //context.strokeStyle = BACKGROUND_COLOR;
  context.globalCompositeOperation = "destination-out";
  context.strokeStyle = "rgba(255,255,255,1)";

  context.lineWidth = lineWidth;
  context.moveTo(x1, y1);
  context.lineTo(x2, y2);
  context.lineCap = "round";
  context.stroke();
  context.closePath();
}



function divideScreen(n) {
  if (n!=3) n=2;
  BoardManager.scrolln=n;
  //console.log("divide the screen ",n)
  let dx=Math.floor(window.innerWidth / n/100)*100;
  var x = container.scrollLeft + dx;
  drawLine(canvas.getContext("2d"), x, 0, x, window.innerHeight, 0.4, BoardManager.getDefaultChalkColor());
  BoardManager.cancelStack.push(BoardManager.getDefaultChalkColor()+","+x+",0,4n"+UI.relative(window.innerHeight,0));
  if (n==3){
    x += dx;
    drawLine(canvas.getContext("2d"), x, 0, x, window.innerHeight, 0.4, BoardManager.getDefaultChalkColor());
    BoardManager.cancelStack.push(BoardManager.getDefaultChalkColor()+","+x+",0,4n"+UI.relative(window.innerHeight,0));
    x += dx;
    drawLine(canvas.getContext("2d"), x, 0, x, window.innerHeight, 0.4, BoardManager.getDefaultChalkColor());
    BoardManager.cancelStack.push(BoardManager.getDefaultChalkColor()+","+x+",0,4n"+UI.relative(window.innerHeight,0));
  }
  BoardManager.saveCurrentScreen();
}


function backgroundClear() {
  canvasBackground.getContext("2d").clearRect(0, 0, canvasBackground.width, canvasBackground.height);
}


function musicScore() {
  backgroundClear();
  let COLORSTAFF = "rgb(128, 128, 255)";
  let fullHeight = window.innerHeight - 32;
  let x = container.scrollLeft;
  let x2 = container.scrollLeft + window.innerWidth;
  let ymiddleScreen = fullHeight / 2;
  let yshift = fullHeight / 7;
  let drawStaff = (ymiddle) => {
    let space = fullHeight / 30;

    for (let i = -2; i <= 2; i++) {
      let y = ymiddle + i * space;
      drawLine(canvasBackground.getContext("2d"), x, y, x2, y, 1.0, COLORSTAFF);
    }
  }


  drawStaff(ymiddleScreen - yshift);
  drawStaff(ymiddleScreen + yshift);

  BoardManager.saveCurrentScreen();
}


let magnetColors = ['rgba(255, 128, 0, 0.9)', 'rgba(0, 128, 0, 0.9)', 'rgba(192, 0, 0, 0.9)', 'rgba(0, 0, 255, 0.9)','rgba(0, 0, 0, 0.9)','rgba(1, 1, 1, 0.1)'];

function nextBackgroundColor(color) {
  for (let i = 0; i < magnetColors.length; i++) {
    if (magnetColors[i] == color) {
      return magnetColors[(i + 1) % magnetColors.length];
    }
  }
  return magnetColors[0];
}


function previousBackgroundColor(color) {
  for (let i = 0; i < magnetColors.length; i++) {
    if (magnetColors[i] == color) {
      return magnetColors[(i - 1) % magnetColors.length];
    }
  }
  return magnetColors[0];
}

/* ***************
   ** MAGNETS   **
   *************** */

class MagnetManager {

  static magnetX = 0;
  static magnetY = 0;
  static currentMagnet = undefined; // last magnet used
  static magnetUnderCursor = undefined;

  static init() {
    document.getElementById("magnetsArrange").onclick = MagnetManager.arrange;
    document.getElementById("magnetsCreateGraph").onclick = MagnetManager.drawGraph;
  }

  static getMagnetUnderCursor() {
    return MagnetManager.magnetUnderCursor;
  }

  /**
   * 
   * @param {boolean} b 
   * @description if b == true, makes all the magnets interactable with the mouse/pointer
   *  if b == false, the magnets cannot be moved
   */
  static setInteractable(b) {
    const v = b ? "auto" : "none";

    let magnets = MagnetManager.getMagnets();

    for (let i = 0; i < magnets.length; i++)
      magnets[i].style.pointerEvents = v;

  }

  /**
   * @returns an array containing all the magnets
   */
  static getMagnets() {
    return document.getElementsByClassName("magnet");
  }


  /**
   * delete all the magnets
   */
  static clearMagnet() {
    MagnetManager.currentMagnet = undefined;
    MagnetManager.magnetX = BoardManager.getCurrentScreenRectangle().x1;
    MagnetManager.magnetY = 0;
    let magnets = MagnetManager.getMagnets();

    while (magnets.length > 0)
      magnets[0].remove();

    Share.sendMagnets();

    Menu.hide();
  }

  /**
   * 
   * @param {*} element 
   * @description add the DOM element element to the list of magnets
   */
  static addMagnet(element) {
    // console.log('addmagnet',element);
    if (MagnetManager.magnetX > BoardManager.getCurrentScreenRectangle().x2 - 10) {
      MagnetManager.magnetX = BoardManager.getCurrentScreenRectangle().x1;
      MagnetManager.magnetY += 64;
    }
    for (;++UI.histcount;){
      if ($id("m"+UI.histcount)==null)
	break;
    }
    //element.id = "m" + Math.random(); //generate randomly an id
    element.id = "m" + UI.histcount;
    ++UI.histcount;
    element.style.left = MagnetManager.magnetX + "px";
    element.style.top = MagnetManager.magnetY + "px";

    //done with setTimeout because images may be loaded
    setTimeout(() => element.style.zIndex = window.innerWidth - element.clientWidth, 400);

    MagnetManager.magnetX += 64;
    MagnetManager.currentMagnet = element;
    element.classList.add("magnet");
    MagnetManager._installMagnet(element);
    //console.log(element);
    document.getElementById("magnets").appendChild(element);
  }


  static arrange() {
    let magnets = MagnetManager.getMagnets();

    for (let i = 0; i < magnets.length; i++) {
      let magnet = magnets[i];
      let x = undefined;
      let y = undefined;

      let magnetContains = (m, x, y) => {
	return (parseInt(m.style.left) <= x && parseInt(m.style.top) <= y &&
		x <= parseInt(m.style.left) + parseInt(m.clientWidth) &&
		y <= parseInt(m.style.top) + parseInt(m.clientHeight));
      }

      let dist = () => {
	let minDist = 100000;
	for (let j = 0; j < magnets.length; j++) {
	  minDist = Math.min(minDist,
			     Math.abs(x - parseInt(magnets[j].style.left)) + Math.abs(y - parseInt(magnets[j].style.top)));
	}
	return minDist;

      }
      let contains = () => {
	for (let j = 0; j < magnets.length; j++) {
	  if (magnetContains(magnets[j], x, y) ||
	      magnetContains(magnets[j], x + parseInt(magnet.clientWidth), y + parseInt(magnet.clientHeight)))
	    return true;
	}
	return false;
      }

      const rect = BoardManager.getCurrentScreenRectangle();

      let generatePosition = () => {
	let count = 0;
	const margin = 32;
	do {
	  x = rect.x1 + (Math.random() * window.innerWidth);
	  y = rect.y1 + (Math.random() * 3 * window.innerHeight / 4);

	  x = Math.max(x, rect.x1 + margin);
	  y = Math.max(y, rect.y1 + margin);
	  x = Math.min(x, rect.x2 - magnet.clientWidth - margin);
	  y = Math.min(y, rect.y2 - magnet.clientHeight - margin);
	  count++;
	}
	while (contains() && count < 50)
      }


      let count = 0;
      let bestDist = 0;
      let bestX = undefined;
      let bestY = undefined;

      while (count < 30) {

	generatePosition();

	if (bestDist < dist()) {
	  bestX = x;
	  bestY = y;
	  bestDist = dist();
	}
	count++;
      }

      magnet.style.left = bestX;
      magnet.style.top = bestY;
    }

  }


  static getNodes() {
    let magnets = MagnetManager.getMagnets();
    let nodes = [];
    for (let i = 0; i < magnets.length; i++) {
      let m = magnets[i];
      nodes.push({ x: parseInt(m.style.left) + parseInt(m.clientWidth) / 2, y: parseInt(m.style.top) + parseInt(m.clientHeight) / 2 });
    }
    console.log(nodes)
    return nodes;
  }

  static drawGraph() {
    MagnetManager.arrange();

    let nodes = MagnetManager.getNodes();
    let context = canvas.getContext("2d");
    let edges = [];
    for (let i = 0; i < nodes.length; i++) {
      edges[i] = [];
      for (let j = 0; j < nodes.length; j++) { edges[i][j] = 0; }
    }

    // returns true iff the line from (a,b)->(c,d) intersects with (p,q)->(r,s)
    function intersects(a, b, c, d, p, q, r, s) {
      var det, gamma, lambda;
      det = (c - a) * (s - q) - (r - p) * (d - b);
      if (det === 0) {
	return false;
      } else {
	lambda = ((s - q) * (r - a) + (p - r) * (s - b)) / det;
	gamma = ((b - d) * (r - a) + (c - a) * (s - b)) / det;
	return (0 < lambda && lambda < 1) && (0 < gamma && gamma < 1);
      }
    };

    let isCrossing = (i, j) => {
      for (let k = 0; k < nodes.length; k++)
	for (let l = 0; l < nodes.length; l++)
	  if (edges[k][l]) {
	    if (intersects(nodes[i].x, nodes[i].y, nodes[j].x, nodes[j].y, nodes[k].x, nodes[k].y, nodes[l].x, nodes[l].y))
	      return true;
	  }
      return false;
    }

    for (let i = 0; i < nodes.length; i++)
      for (let j = 0; j < nodes.length; j++) {
	if (Math.abs(nodes[i].x - nodes[j].x) + Math.abs(nodes[i].y - nodes[j].y) < 400 && !isCrossing(i, j)) {
	  edges[i][j] = 1;
	  drawLine(context, nodes[i].x, nodes[i].y, nodes[j].x, nodes[j].y);
	}
      }

    BoardManager.save();
  }


  static _installMagnetsNoMsg() {

    let magnets = MagnetManager.getMagnets();

    for (let i = 0; i < magnets.length; i++)
      MagnetManager._installMagnet(magnets[i]);

  }


  static installMagnets() {
    MagnetManager._installMagnetsNoMsg();
    Share.sendMagnets();
  }

  static _installMagnet(element) {
    makeDraggableElement(element);
    var divText=element.firstChild;
    if (divText!=null){
      divText.onpointerdown = (e) => {  e.stopPropagation(); }
      divText.onpointermove = (e) => { e.stopPropagation(); }
      divText.onpointerup = (e) => {  e.stopPropagation(); }
      divText.onkeydown = (e) => { UI.textmagnet_onkeydown(e,divText); }
    }
    element.onmouseenter = () => { MagnetManager.magnetUnderCursor = element };
    element.onmouseleave = () => { MagnetManager.magnetUnderCursor = undefined };

    function makeDraggableElement(element) {
      //console.log('makedraggable',element.id);
      var pos1 = 0, pos2 = 0, pos3 = 0, pos4 = 0;
      
      element.addEventListener("pointerdown", dragMouseDown);

      TouchScreen.addTouchEvents(element);

      let otherElementsToMove = [];
      let canvasCursorStore = undefined;
      let drag = true;

      function dragMouseDown(e) {
	//console.log('dragmousedown',element);
	if (UI.focused!=0)
	  return;
	drag = true;
	/**
	 * 
	 * @param {*} el 
	 * @param {*} bigel 
	 * @returns true if el is inside bigel
	 */
	function inside(el, bigel) {
	  return el.offsetLeft > bigel.offsetLeft && el.offsetTop > bigel.offsetTop &&
	    el.offsetLeft + el.clientWidth < bigel.offsetLeft + bigel.clientWidth &&
	    el.offsetTop + el.clientHeight < bigel.offsetTop + bigel.clientHeight;
	}
	var canvas=$id("canvas");
	canvasCursorStore = canvas.style.cursor;
	e = e || window.event;
	e.preventDefault(); //to avoid the drag/drop by the browser
	// get the mouse cursor position at startup:
	pos3 = e.clientX;
	pos4 = e.clientY;
	document.onpointerup = closeDragElement;
	document.onmouseup = closeDragElement;
	// call a function whenever the cursor moves:
	document.onpointermove = elementDrag;

	let magnets = MagnetManager.getMagnets();
	otherElementsToMove = [];

	//if(elmt.style.clipPath == undefined) //if not an image (otherwise bug)
	for (let i = 0; i < magnets.length; i++)
	  if (magnets[i] != element && inside(magnets[i], element)) {
	    otherElementsToMove.push(magnets[i]);
	  }


      }




      function elementDrag(e) {
	if (!drag) return;
	if (e.target.parentNode && e.target.parentNode.classList.contains("magnetText")){
	  MagnetManager.currentMagnet = e.target.parentNode;
	  e.target.parentNode.classList.add("magnetDrag");
	}
	else {
	  MagnetManager.currentMagnet = e.target;
	  if (e.target.classList) e.target.classList.add("magnetDrag");
	}

	canvas.style.cursor = "none";
	e = e || window.event;
	e.preventDefault();
	// calculate the new cursor position:
	pos1 = pos3 - e.clientX;
	pos2 = pos4 - e.clientY;
	pos3 = e.clientX;
	pos4 = e.clientY;



	// set the element's new position:
	Share.execute("magnetMove", [element.id, element.offsetLeft - pos1, element.offsetTop - pos2]);

	for (let el of otherElementsToMove) {
	  Share.execute("magnetMove", [el.id, el.offsetLeft - pos1, el.offsetTop - pos2]);
	}
      }

      function closeDragElement(e) {
	if (!drag)
	  return;
	drag = false;
	if (!e.target.parentNode)
	  return;
	console.log("close drag",e.target.parentNode.style.top);// alert("close drag");
	UI.link(0);
	e.target.classList.remove("magnetDrag");
	canvas.style.cursor = canvasCursorStore;

	// stop moving when mouse button is released:
	document.onmouseup = null;
	document.onmousemove = null;
	let f=e.target;
	if (!f.classList.contains("magnet"))
	  f=f.parentNode;
	//console.log(f);
	if (parseInt(f.style.top)<(Menu.isShown()?10:0)
	   ){
	  f.style.visibility='hidden';
	  UI.ck_display_trash_msg();
	  //console.log('hide',f,f.style.visibility);
	  //f.remove();
	}
      }
    }
  }




  static addMagnetImage(filename,x=0,y=0) {
    let content='<img src="magnets/'+filename+'">'
    UI.addmagnet(x,y,'comment',content,-1);
    return;
    let img = new Image();
    img.src = "magnets/" + filename;
    MagnetManager.addMagnet(img);
  }


  static addMagnetText(x, y,txt="textarea",evaluator='',content='',evalmode=0) {
    //console.log('addmagnettext',x,y,evaluator,UI.magnet_evaluator,container.scrollLeft,window.innerWidth);
    if (evaluator==''){
      evaluator=UI.magnet_evaluator;
    }
    let div = document.createElement("div");
    MagnetManager.addMagnet(div);
    if (x-container.scrollLeft>window.innerWidth-300){
      //console.log('shift new magnet position',x,container.scrollLeft+window.innerWidth-300);
      x=container.scrollLeft+window.innerWidth-300;
    }
    div.style.left = (x-container.scrollLeft) + "px";
    div.style.top = y + "px";
    div.classList.add("magnetText");
    var col= palette.colors[0],bg=col=="white"?"rgba(0.2, 0.2, 0.2, 0.9)":"white";
    //console.log('Palette col bg',col,bg);
    div.style.backgroundColor = bg;
    div.style.color = col;
    if (evaluator=='svg' ){
      //console.log(evaluator,content);
      div.innerHTML=content;
      //console.log(div.firstChild.clientWidth);
      div.style.width=div.firstChild.clientWidth;
      div.style.height=div.firstChild.clientHeight;
      return div;
    }
    if (evaluator=='img'){
      let img = new Image();
      img.src = content;
      div.appendChild(img);
      return div;
    }
    if (evaluator=='comment'){
      //console.log('addmagnettext',x,y);
      let s=UI.addcomment(content);
      div.innerHTML=s;
      div.classList.add(evaluator);
      if (evalmode==0){
	div.firstChild.nextSibling.click();
	UI.tableaucm=1;
      }
      UI.link();
      return div;
    }
    let divtype = document.createElement("span");
    divtype.style.fontSize=UI.detectmob()?'16px':'20px';
    div.appendChild(divtype);
    if (evaluator=='handwriting'){
      //console.log('handw',content);
      divtype.classList.add(evaluator);
      divtype.innerText=content;
      divtype.style.display='none';
      let c=UI.handwrite(content);
      div.appendChild(c);
      UI.link();
      return div;
    }
    if (evaluator=='xcas' || evaluator=='cas'){
      UI.ckswitch(evaluator,0);
      divtype.innerHTML='>';
    }
    if (evaluator=='python' || evaluator=='py'){
      UI.ckswitch(evaluator,0);
      divtype.innerHTML='>>>';
    }
    if (evaluator=='javascript' || evaluator=='js'){
      UI.ckswitch(evaluator,0);      
      divtype.innerHTML='>>';
    }
    let divText = document.createElement((txt=="cm")?'textarea':txt);
    //console.log('classlist.add',evaluator); console.trace();
    divText.classList.add(evaluator);
    div.appendChild(divText);
    if (UI.detectmob()){
      divText.style.fontSize=16+'px';
      div.style.fontSize=16+'px';
    }
    else {
      divText.style.fontSize=20+'px';
      div.style.fontSize=20+'px';
    }
    // let divText = document.createElement("div");
    if (txt=="cm" || txt=="textarea" || txt=="textinput"){
       if (content=='') {
	 if (UI.micropy==-1)
	   divText.innerHTML="JS";
	 if (UI.micropy==0)
	   divText.innerHTML="CAS";
	 if (UI.micropy==1)
	   divText.innerHTML="Py";
	 divText.innerHTML=" ";
       }
       else 
	 divText.innerHTML=content;
    }
    if (txt=="div")
      divText.innerHTML = "x^2";
    if (1 || txt!="cm"){
      divText.onpointerdown = (e) => { e.stopPropagation(); }
      divText.onpointermove = (e) => { e.stopPropagation(); }
      divText.onpointerup = (e) => { e.stopPropagation(); }
      divText.onkeydown = (e) => { UI.textmagnet_onkeydown(e,divText); }
    }
    if (!UI.detectmob() && x>canvas.width-512)
      x=canvas.width-512; // should be something like parseInt(div.style.maxWidth);
    divText.style.backgroundColor = bg;
    divText.style.color = col;
    divText.contentEditable = true;
    divText.style.height=(parseInt(divText.style.fontSize)+8)+"px";
    if (evalmode==0){
      divText.focus();
      //document.execCommand('selectAll', false, null);
    }
    if (txt=="textarea"){
      //divText.select();
      return divText;
    }
    if (txt=="cm"){
      //divText.select();
      divText.addEventListener('focus',()=>{UI.textareatocm(divText);});
      if (evalmode==0){
	let cm=UI.textareatocm(divText);
	// cm.execCommand('selectAll');
	//console.log(cm);
      }
      if (evalmode==1){
	//console.log('evalmode',divText);
	UI.tableaueval(divText);
	//console.log(divText.type);
      }
      if (evalmode==1 || evalmode==-1){
	UI.resizetextarea(divText); let h=(2+parseInt(divText.style.fontSize))*divText.rows+'px';divText.style.height=h; // console.log('addmagnettext',divText.rows,h,divText.style.height);
      }
    }
    return divText;
  }



  static removeCurrentMagnet() {
    if (MagnetManager.currentMagnet == undefined)
      return;
    MagnetManager.currentMagnet.style.visibility='hidden';//remove();
    MagnetManager.currentMagnet == undefined;
    MagnetManager.magnetUnderCursor = undefined;
  }

  static printCurrentMagnet() {
    const img = MagnetManager.currentMagnet;

    if (!(img instanceof Image)) {
      // fixme: handwriting
      let cur=MagnetManager.currentMagnet;
      //console.log(cur);
      if (cur.tagName=='DIV' && cur.hasChildNodes()){
	let mx=parseInt(cur.style.left)+container.scrollLeft,my=parseInt(cur.style.top);
	if (cur.firstChild.nextSibling && cur.firstChild.nextSibling.tagName=='CANVAS'){
	  let s_in = cur.firstChild.innerText;
	  // console.log('magnet hand',mx,my,s_in);
	  if (s_in.length && s_in[s_in.length-1]=='&')
	    s_in=s_in.substr(0,s_in.length-1);
	  let tab=s_in.split('&');
	  let dims=UI.find_whxy(tab);
	  mx += dims[2]; my += dims[3];
	  // console.log('magnet hand',tab,mx,my);
	  let ctx=canvas.getContext("2d");
	  for (let i=0;i<tab.length;++i){
	    // tab[i] will be translated by mx,my and pushed in cancelstack
	    let lxypcs=UI.decode_polygon(tab[i]);
	    let s=lxypcs[3]=='white'?'':(lxypcs[3]+',');
	    s += (lxypcs[0][0]+mx)+','+(lxypcs[1][0]+my)+','+lxypcs[4];
	    BoardManager.cancelStack.push(s);
	    // console.log(tab[i],mx,my,s,BoardManager.cancelStack.stack.length,BoardManager.cancelStack.currentIndex,BoardManager.cancelStack.stack);
	  }
	  BoardManager.redraw_handwritings();
	  return;
	}
	let s=cur.innerHTML;
	let copy=document.createElement('DIV');
	copy.innerHTML=s;
	copy.style.backgroundColor=cur.style.backgroundColor;
	copy.style.color=cur.style.color;
	MagnetManager.addMagnet(copy);
	return;
      }
      MagnetManager.addMagnet(cur);
      console.log("the current magnet could not be printed!")
      return;
    }

    let context = document.getElementById("canvas").getContext("2d");

    let x = parseInt(img.style.left);
    let y = parseInt(img.style.top);
    let s = img.style.clipPath;

    s = s.substr("polygon(".length, s.length - "polygon(".length - ")".length);

    context.globalCompositeOperation = "source-over";
    context.save();
    context.beginPath();
    let begin = true;
    for (let pointStr of s.split(",")) {
      pointStr = pointStr.trim();
      let a = pointStr.split(" ");
      if (begin)
	context.moveTo(x + parseInt(a[0]), y + parseInt(a[1]));
      else
	context.lineTo(x + parseInt(a[0]), y + parseInt(a[1]));
      begin = false;
    }
    context.closePath();
    context.clip();

    context.drawImage(img, x, y);

    context.restore();



    BoardManager.save();
  }


  static register(magnetSetName) {
    document.getElementById(magnetSetName).onclick = eval(magnetSetName);
  }


}


/** Gale-Shapley */

function createMagnet(content) {
  let o = document.createElement("div");
  o.innerHTML = content;
  return o;
}


function magnetsClear() {
  Share.execute("magnetsClear", []);
}


function magnetGS() {
  Welcome.hide(); Menu.hide();
  for (let i=1;i<=3;++i){
    UI.addmagnet(40*i,0,'comment',''+i,-1);
    let o=MagnetManager.currentMagnet;
    o.style.backgroundColor='white';
    o.style.color='black';
  }
  for (let i=1;i<=3;++i){
    UI.addmagnet(120+40*i,0,'comment',''+i,-1);
    let o=MagnetManager.currentMagnet;
    o.style.backgroundColor='white';
    o.style.color='orange';
    o.classList.add("GS_B")
  }
  UI.link();
  return;

}

/** Sorting */

function createMagnetRainbow(i) {
  UI.addmagnet(40*i,0,'comment',''+i,-1);
  let o = MagnetManager.currentMagnet;//createMagnet(content);
  let colors = ['rgb(139, 97, 195)', 'rgb(115, 97, 195)', 'rgb(93, 105, 214)', 'rgb(40, 167, 226)', 'rgb(40, 204, 226)', 'rgb(40, 226, 201)', 'rgb(40, 226, 148)',
		'rgb(40, 226, 102)', 'rgb(130, 226, 40)', 'rgb(170, 226, 40)', 'rgb(223, 226, 40)', 'rgb(226, 183, 40)',
		'rgb(226, 152, 40)', 'rgb(226, 124, 40)', 'rgb(226, 77, 40)', 'rgb(255, 0, 0)', 'rgb(144, 24, 24)'];
  o.style.backgroundColor = colors[i - 1];
  return o;
}


function magnetSorting() {
  Welcome.hide(); Menu.hide();
  for (let i = 1; i <= 17; i++)
    createMagnetRainbow(i)
  UI.link();
  UI.set_readonly(true);
}

/** B-trees */

function magnetBTrees() {
  Welcome.hide(); Menu.hide();
  for (let i = 1; i <= 17; i++)
    createMagnetRainbow(i);

  for (let i = 0; i < 7; i++)
    MagnetManager.addMagnetImage("Btreenode.png",50+Math.floor(i/3)*300,50+(i%3)*150);

  UI.set_readonly(true);

}

/** Graphs */

function magnetGraphNodes() {
  Welcome.hide(); Menu.hide();
  //MagnetManager.clearMagnet();
  for (let i of ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']){
    UI.addmagnet(40*(i.charCodeAt(0)-64),0,'comment',i,-1);
    let o = MagnetManager.currentMagnet;
    o.style.backgroundColor='white';
    o.style.color='black';    
  }
  UI.set_readonly(true);
}

/** magnetFloydsAlgorithm **/
function magnetFloydsAlgorithm() {
  Welcome.hide(); Menu.hide();
  MagnetManager.addMagnetImage("turtlerabbit/turtle.png",50,0);
  MagnetManager.addMagnetImage("turtlerabbit/rabbit.png",200,0);
  
}

function magnetGraphSimCity() {
  Welcome.hide(); Menu.hide();
  //MagnetManager.clearMagnet();
  let i=0;
  const simCityPictures = ["antenne.png",       "commerce.png",      "parking.png",        "tour.png",
			   "batimentplat.png",  "foursolaire.png",   "residence2.png",    "usine.png",
			   "building.png",      "gare.png",          "residencebleu.png",
			   "chateaudeau.png",   "nuclearplant.png",  "residence.png",
			   "citerne.png",       "parc.png",          "stade.png"
			  ];

  for(let name of simCityPictures) {
    MagnetManager.addMagnetImage("simCityGraph/" + name,50+Math.floor(i/3)*300,50+(i%3)*150);
    ++i;
  }
}

/** Tilings */

function createTiling(leftColor, upColor, rightColor, bottomColor, leftText, upText, rightText, bottomText ) {
  let xmlns = "http://www.w3.org/2000/svg";
  var div = document.createElement("div");
  let size = 100;
  var svgElem = document.createElementNS(xmlns, "svg");
  svgElem.setAttributeNS(null, "viewBox", "0 0 " + size + " " + size);
  svgElem.setAttributeNS(null, "width", size);
  svgElem.setAttributeNS(null, "height", size);
  svgElem.style.display = "block";
  

  function createPath(pathDesc, color) {
    var path = document.createElementNS(xmlns, "path");
    path.setAttributeNS(null, 'stroke', "#333333");
    path.setAttributeNS(null, 'stroke-width', 10);
    path.setAttributeNS(null, 'stroke-linejoin', "round");
    path.setAttributeNS(null, 'd', pathDesc);
    path.setAttributeNS(null, 'fill', color);
    path.setAttributeNS(null, 'opacity', 1.0);
    return path;
  }
  
  svgElem.appendChild(createPath("M 50 50 L 0 0 L 0 100 Z", leftColor));
  svgElem.appendChild(createPath("M 50 50 L 0 0 L 100 0 Z", upColor));
  svgElem.appendChild(createPath("M 50 50 L 100 0 L 100 100 Z", rightColor));
  svgElem.appendChild(createPath("M 50, 50 L 100 100 L 0 100 Z", bottomColor));
  div.appendChild(svgElem);
  div.style.padding = "0px";
  div.style.width=size+'px';
  div.style.height=size+'px';
  return div;

}


function magnetTilings() {
  Welcome.hide(); Menu.hide();
  MagnetManager.addMagnet(createTiling("yellow", "red", "green", "red"));
  MagnetManager.addMagnet(createTiling("green", "red", "green", "yellow"));
  MagnetManager.addMagnet(createTiling("green", "red", "green", "yellow"));

  MagnetManager.addMagnet(createTiling("red", "red", "red", "red"));
  MagnetManager.addMagnet(createTiling("red", "yellow", "red", "green"));
  MagnetManager.addMagnet(createTiling("red", "yellow", "yellow", "yellow"));

  MagnetManager.addMagnet(createTiling("green", "red", "green", "yellow"));
  MagnetManager.addMagnet(createTiling("green", "green", "red", "green"));
  MagnetManager.addMagnet(createTiling("red", "yellow", "red", "green"));
  UI.link();
}



/**
 *  picture not available
 */
function magnetUnionFind() {
  Welcome.hide(); Menu.hide();
  MagnetManager.addMagnetImage("unionfind0.png");
}

function loadMagnets() {
  MagnetManager.register("magnetGS");
  MagnetManager.register("magnetSorting");
  MagnetManager.register("magnetBTrees");
  MagnetManager.register("magnetGraphNodes");
  MagnetManager.register("magnetTilings");
  // MagnetManager.register("magnetUnionFind");
  MagnetManager.register("magnetGraphSimCity");
  MagnetManager.register("magnetFloydsAlgorithm");
}


/*  ****************
    ** UNDO STACK **
    ****************  */

class CancelStack {
  stack = [];
  currentIndex = 0;
  n = 0;

  clear() {
    this.stack = [];
    this.currentIndex = 0;
    this.n = 0;
  }

  push(data) {
    this.stack[this.currentIndex] = data;
    this.currentIndex++;
    this.n = this.currentIndex;
    //console.log('cancelstack push',data,this.n);console.trace();
  }


  back(steps=1) {
    if (steps>0 && this.currentIndex>=steps)
      this.currentIndex -= steps;
    else
      this.currentIndex = 0;
    //console.log('back',this.stack.slice(0,this.currentIndex));
    return this.stack.slice(0,this.currentIndex);
  }

  forward(steps=1) {
    if (steps>0 && this.currentIndex+steps<=this.n)
      this.currentIndex += steps;
    else 
      this.currentIndex = this.n;     
    return this.stack.slice(0,this.currentIndex);
  }

}

/**
  * Manage the board
  */
class BoardManager {

  /** name of the board. Default is 0 (this name is used for storing in localStorage) */
  static boardName = 0;
  static scrolln=3;

  /** stack to store the cancel/redo actions */
  static cancelStack = new CancelStack();

  static _rightExtendCanvasEnable = true;

  /**
   * initialization (button)
   */
  static init() {
    document.getElementById("blackboardClear").onclick = () => {
      BoardManager._clear();
      BoardManager.save();
      Menu.hide();
    }
  }


  static getBackgroundColor() {
    return document.getElementById("canvasBackground").style.backgroundColor;
  }
  
  /**
   * erase the board
   */
  static _clear() {
    document.getElementById("canvas").width = document.getElementById("canvas").width; //clear
    BoardManager.cancelStack.clear();
  }


  /**
   * resize the board to the size of the window
   */
  static resize() {
    if (document.getElementById("canvas").width < window.innerWidth)
      document.getElementById("canvas").width = window.innerWidth;

    if (document.getElementById("canvas").height < window.innerHeight)
      document.getElementById("canvas").height = window.innerHeight;
  }


  static getDefaultChalkColor() {
    return (document.getElementById("canvasBackground").style.backgroundColor == "black") ? "white" : "black";
  }


  static _createCanvasForRectangle(r) {
    let C = document.createElement("canvas");
    C.width = r.x2 - r.x1;
    C.height = r.y2 - r.y1;
    let ctx = C.getContext("2d");
    ctx.drawImage(document.getElementById("canvas"),
		  r.x1, r.y1, r.x2 - r.x1, r.y2 - r.y1, //coordinates in the canvas
		  0, 0, r.x2 - r.x1, r.y2 - r.y1); //coordinates in the magnet
    return C;
  }


  /**
   * 
   * @param {*} r a rectangle {x1, y1, x2, y2}
   * @returns the content as a string of the image
   */
  static _toBlobOfRectangle(r, callback) {
    return BoardManager._createCanvasForRectangle(r).toBlob(callback);
  }



  /** disabled
   * save the current board into the cancel/redo stack but also in the localStorage of the browser
   */
  static save(rectangle) {
    return; 
    // if (rectangle == undefined) {
    document.getElementById("canvas").toBlob((blob) => {
      console.log("save that blob: " + blob)
      //  localStorage.setItem(Share.getTableauNoirID(), canvas.toDataURL());
      BoardManager.cancelStack.push(blob);
      //Share.sendFullCanvas(blob);
    }
					    );
    /*}
      else {
      BoardManager._toBlobOfRectangle(rectangle, (blob) => {
      rectangle.blob = blob;
      BoardManager.cancelStack.push(rectangle);
      });
      
      document.getElementById("canvas").toBlob((blob) => {
      console.log("save that blob: " + blob)
      localStorage.setItem(BoardManager.boardName, blob);
      }
      );
      }
    */


  }




  static getCurrentScreenRectangle() {
    const x1 = container.scrollLeft;
    const y1 = container.scrollTop;
    const x2 = x1 + window.innerWidth;
    const y2 = y1 + window.innerHeight;
    return { x1: x1, y1: y1, x2: x2, y2: y2 };
  }

  static saveCurrentScreen() {
    BoardManager.save(BoardManager.getCurrentScreenRectangle());
  }

  /**
   * load the board from the local storage
   */
  static load(data = localStorage.getItem(Share.getTableauNoirID())) {
    // let data = localStorage.getItem(BoardManager.boardName);

    if (data != undefined) {
      BoardManager._clear();
      try {
        let image = new Image();
        image.src = data;
        image.onload = function () {
          document.getElementById("canvas").width = image.width;
          document.getElementById("canvas").height = image.height;
          document.getElementById("canvas").getContext("2d").drawImage(image, 0, 0);
          BoardManager.save();
          console.log("loaded!")
        }
      }
      catch (e) {

      }

    }
    else {
      BoardManager._clear();
      BoardManager.save();
    }

  }




  /**
   * load the board from the local storage
   */
  static loadWithoutSave(data = localStorage.getItem(BoardManager.boardName)) {
    // let data = localStorage.getItem(BoardManager.boardName);

    if (data != undefined) {
      BoardManager._clear();
      let image = new Image();
      image.src = data;
      image.onload = function () {
        document.getElementById("canvas").width = image.width;
        document.getElementById("canvas").height = image.height;
        document.getElementById("canvas").getContext("2d").drawImage(image, 0, 0);
        console.log("loaded!")
      }
    }
    else {
      BoardManager._clear();
    }

  }

  /**
   * @returns the number of pixels when scrolling
   */
  static scrollQuantity() {
    let dx= Math.floor(window.innerWidth / BoardManager.scrolln/100)*100;
    if (UI.detectmob() || dx<=200)
      return 200;
    return dx;
  }

  static shift_magnets(dx){
    let magnets = document.getElementsByClassName("magnet");
    for (let i = 0; i < magnets.length; i++) {
      // console.log('right', magnets[i].style.left,BoardManager.scrollQuantity());
      var tmp=parseInt(magnets[i].style.left);
      tmp += dx;
      magnets[i].style.left = tmp+"px";
    }
  }

  /**
   * go left
   */
  static left() {
    if (UI.locked){ console.log('locked'); return ; }
    // BoardManager.cancelStack.clear();     BoardManager.save();
    let dx=BoardManager.scrollQuantity();
    let x = container.scrollLeft - dx;
    if (x<0 && x>-1)
      x=0;
    if (x < 0) {
      BoardManager.showPageNumber(0);
      return;
    }

    UI.locked=true;
    container.scrollTo({ top: 0, left: x }); // behavior: 'smooth'
    let magnets = document.getElementsByClassName("magnet");
    //console.log('left',dx);
    for (let i = 0; i < magnets.length; i++) {
      // console.log('left', magnets[i].style.left,BoardManager.scrollQuantity());
      var tmp=parseInt(magnets[i].style.left);
      tmp += dx;
      magnets[i].style.left = tmp+"px";
    }

    user.x -= dx;
    BoardManager.showPageNumber(x);
    UI.locked=false;
  }

  /**
   * go right (and extend the board if necessary)
   */
  static right() {
    if (UI.locked){ console.log('locked'); return ; }
    // BoardManager.cancelStack.clear();     BoardManager.save();
    UI.tableaulocalsave('tableaunoirxcas');
    const MAXCANVASWIDTH = 16384;

    if (container.scrollLeft >= MAXCANVASWIDTH - window.innerWidth) {
      console.log('maxscanwidth');
      container.scrollLeft = MAXCANVASWIDTH - window.innerWidth;
      return;
    }
    UI.locked=true;
    if ((container.scrollLeft >= canvas.width - window.innerWidth - BoardManager.scrollQuantity()) && BoardManager._rightExtendCanvasEnable) {
      let image = new Image();
      image.src = canvas.toDataURL();
      console.log("extension: canvas width " + canvas.width + " to " + (container.scrollLeft + window.innerWidth))
      canvas.width = ((canvas.width / BoardManager.scrollQuantity()) + 1) * BoardManager.scrollQuantity();
      BoardManager.redraw_handwritings();
      /*
      const context = document.getElementById("canvas").getContext("2d");
      context.globalCompositeOperation = "source-over";
      context.globalAlpha = 1.0;
      image.onload = function () {
        context.drawImage(image, 0, 0);
      }
      BoardManager._rightExtendCanvasEnable = false;
      setTimeout(() => { BoardManager._rightExtendCanvasEnable = true }, 1000);
      */
    }
    const x = container.scrollLeft + BoardManager.scrollQuantity();
    console.log('right',container.scrollLeft,x);
    UI.scroll_magnets(BoardManager.scrollQuantity());
    container.scrollTo({ top: 0, left: x }); // , behavior: 'smooth'
    user.x += BoardManager.scrollQuantity();
    BoardManager.showPageNumber(x);
    UI.locked=false;
  }


  static showPageNumber(x) {

    pageNumber.classList.remove("pageNumberHidden");
    pageNumber.classList.remove("pageNumber");
    setTimeout(() => {
      const n = Math.round(x / BoardManager.scrollQuantity());
      const total = Math.round(canvas.width / BoardManager.scrollQuantity());
      container.scrollLeft = (n) * BoardManager.scrollQuantity();
      pageNumber.innerHTML = (n + 1) + "/" + (total); pageNumber.classList.add("pageNumber");
    }, 300)

  }

  /**
   * 
   * @param {*} level 
   */
  static _loadCurrentCancellationStackData(data) {
    // console.log(data.length);
    let ctx = document.getElementById("canvas").getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (let i=0;i<data.length;++i){
      if (data[i]!=undefined)
	UI.draw_encoded(data[i],ctx);
    }
    return;
    /* old code */
    let image = new Image();

    const context = document.getElementById("canvas").getContext("2d");
    context.globalCompositeOperation = "source-over";
    context.globalAlpha = 1.0;

    //  if (data instanceof Blob) {
    image.src = URL.createObjectURL(data);
    image.onload = function () {
      document.getElementById("canvas").width = image.width;
      document.getElementById("canvas").height = image.height;
      context.drawImage(image, 0, 0);
    }
    /*  }
        else {
        console.log("_loadCurrentCancellationStackData with rectangle at " + data.x1)
        image.src = URL.createObjectURL(data.blob);
        image.onload = function () {
        context.clearRect(data.x1, data.y1, data.x2 - data.x1, data.y2 - data.y1);
        context.drawImage(image, data.x1, data.y1);
        }
        }*/ //still bugy


  }

  static redraw_handwritings(){
    //console.log('redraw_handwritings',BoardManager.cancelStack.stack.length,BoardManager.cancelStack.currentIndex);
    let data=BoardManager.cancelStack.stack.slice(0,BoardManager.cancelStack.currentIndex);
    //console.log('redraw_handwritings',''+data,BoardManager.cancelStack.currentIndex);
    this._loadCurrentCancellationStackData(data);
  }

  /**
   * 
   */
  static cancel(steps=1) {
    if (!Share.isShared())
      BoardManager._loadCurrentCancellationStackData(BoardManager.cancelStack.back(steps));
  }



  /**
   * 
   */
  static redo(steps=1) {
    if (!Share.isShared())
      BoardManager._loadCurrentCancellationStackData(BoardManager.cancelStack.forward(steps));
  }
}

class ChalkCursor {

  /**
   * 
   * @param {*} color 
   * @returns the .style.cursor of the canvas if you want to have a cursor that looks like a chalk with the color color.
   */
  static getStyleCursor(color) {
    return `url(${ChalkCursor.getCursorURL(color)}) 0 0, auto`;
  }


  static getCursorURL(color) {
    const sizeX = 32;
    const sizeY = 44;
    const chalkRadius = 3;
    const canvas = document.createElement('canvas');
    canvas.width = sizeX;
    canvas.height = sizeY;
    const context = canvas.getContext("2d");

    /*  ctx.beginPath();
        ctx.moveTo(0, 0);
        ctx.lineTo(sizeX, sizeY);
        ctx.lineWidth = chalkRadius * 2;
        ctx.strokeStyle = color;
        ctx.stroke();
	
        ctx.beginPath();
        ctx.arc(chalkRadius, chalkRadius, chalkRadius, 0, 2 * Math.PI);
        ctx.strokeStyle = "black";
        ctx.lineWidth = 1;
        ctx.stroke();*/

    const angle = Math.atan2(sizeY, sizeX);
    const angleOpening = 0.3;
    const anglePlus = angle + angleOpening;
    const angleMinus = angle - angleOpening;
    const sizeHead = 16;
    const length = 34;
    const p1 = { x: sizeHead * Math.cos(anglePlus), y: sizeHead * Math.sin(anglePlus) };
    const p2 = { x: sizeHead * Math.cos(angleMinus), y: sizeHead * Math.sin(angleMinus) };
    const ll = {x : length * Math.cos(angle), y: length * Math.sin(angle) };

    context.beginPath();
    context.moveTo(0, 0);
    context.lineTo(p1.x, p1.y);
    context.lineTo(p1.x + ll.x, p1.y + ll.y);
    context.lineTo(p2.x + ll.x, p2.y + ll.y);
    context.lineTo(p2.x, p2.y);
    context.lineTo(0, 0);

    context.lineWidth = 1;
    context.strokeStyle = "black";
    context.stroke();
    context.fillStyle = color;
    context.fill();

    context.beginPath();
    context.moveTo(sizeHead * Math.cos(anglePlus), sizeHead * Math.sin(anglePlus));
    context.lineTo(sizeHead * Math.cos(angleMinus), sizeHead * Math.sin(angleMinus));
    context.stroke();

    return canvas.toDataURL();
  }
}

class EraserCursor {

  // `url('img/eraser.png') 0 0, auto`;

  /**
   * 
   * @param {*} size 
   * @returns the .style.cursor of the canvas if you want to have a cursor that looks like an eraser of size size
   */
  static getStyleCursor(size = 20) {
    if(size > 128) size = 128;
    return `url(${EraserCursor.getCursorURL(size)}) ${size/2} ${size/2}, auto`;
  }


  static getCursorURL(size) {
    const radius = size / 2;
    const canvas = document.createElement('canvas');
    canvas.width = 2 * radius;
    canvas.height = 2 * radius;
    const ctx = canvas.getContext("2d");

    ctx.beginPath();
    ctx.arc(radius, radius, radius, 0, 2 * Math.PI);
    ctx.strokeStyle = BoardManager.getBackgroundColor() == "black" ? "white" : "black";
    ctx.lineWidth = 1;
    ctx.stroke();
    ctx.fillStyle =  BoardManager.getBackgroundColor() == "black" ? "rgba(255, 255, 255, 0.2)" : "rgba(0, 0, 0, 0.2)";
    ctx.fill();
    return canvas.toDataURL();
  }
}

class LoadSave {

  /**
   * @description initialize the button Save and the drag and drop of files
   */
  static init() {
    //document.getElementById("save").onclick = LoadSave.save;
    //document.getElementById("load").onchange = LoadSave.load;
    document.body.ondragover = (event) => {
      // Prevent default behavior (Prevent file from being opened)
      event.preventDefault();
    }
    document.body.ondrop = (event) => {
      // Prevent default behavior (Prevent file from being opened)
      event.preventDefault();

      if (event.dataTransfer.items) {
        // Use DataTransferItemList interface to access the file(s)
        for (var i = 0; i < event.dataTransfer.items.length; i++) {
          // If dropped items aren't files, reject them
          if (event.dataTransfer.items[i].kind === 'file') {
            var file = event.dataTransfer.items[i].getAsFile();
            LoadSave.loadFile(file);
          }
        }
      } else {
        // Use DataTransfer interface to access the file(s)
        for (var i = 0; i < event.dataTransfer.files.length; i++) {
          LoadSave.loadFile(file[i]);
        }
      }
      console.log('drop');
      Menu.hide(); Welcome.hide();
      window.setTimeout(()=>{UI.link();},2000);

    };
  }

  /**
   * 
   * @param {File} file 
   * @description load the file file
   */
  static loadFile(file) {
    if (file) {
      let reader = new FileReader();
      reader.onerror = function (evt) { }

      /** load a .tableaunoir file, that is, a file containing the blackboard + some magnets */
      if (file.name.endsWith(".tableaunoir")) {
	container.scrollTo({ top: 0, left: 0 }); // behavior: 'smooth'
        reader.readAsText(file, "UTF-8");
        reader.onload = function (evt) { LoadSave.loadJSON(JSON.parse(evt.target.result)); }
      }
      else {
        /** load an image and add it as a magnet */
        reader.readAsDataURL(file);
        reader.onload = function (evt) {
          let img = new Image();
          img.src = evt.target.result;
          MagnetManager.addMagnet(img);
        }
      }
    }
  }



  /**
   * 
   * @param {*} obj 
   * @description load the JSON object:
   * obj.canvasDataURL is the content of the canvas
   * obj.magnets is the HTML code of the magnets
   */
  static loadJSON(obj) {
    BoardManager.load(obj.canvasDataURL);
    console.log('load magnets',obj.magnets);
    document.getElementById("magnets").innerHTML = obj.magnets;
    MagnetManager.installMagnets();
    let magnets = document.getElementsByClassName("magnet");
    UI.tableaufixload(magnets);
  }

  /**
   * @description save the blackboard and the magnets
   */
  static save() {
    const x1 = container.scrollLeft;
    console.log('scrollleft',x1);
    BoardManager.shift_magnets(x1); // move magnets to their position if container not scrolled
    let magnets = document.getElementById("magnets").innerHTML;
    BoardManager.shift_magnets(-x1);
    console.log('save magnets',magnets);
    let canvasDataURL = document.getElementById("canvas").toDataURL();
    let obj = { magnets: magnets, canvasDataURL: canvasDataURL };
    let js=JSON.stringify(obj);
    if (window.localStorage){
      console.log('saved in localStorage',document.getElementById("name").value);
      localStorage.setItem(document.getElementById("name").value,js);
    }
    LoadSave.download(document.getElementById("name").value + ".tableaunoir", js);
  }

  /**
   * 
   * @param {*} filename 
   * @param {*} text 
   * @description propose to download a file called filename that contains the text text
   */
  static download(filename, text) {
    let element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
  }

}

class Menu {
  static init() {
    // Get the element with id="defaultOpen" and click on it
    document.getElementById("defaultOpen").click();
  }

  static openPage(pageName, elmnt) {
    // Hide all elements with class="tabcontent" by default */
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
      tabcontent[i].style.display = "none";
    }

    // Remove the background color of all tablinks/buttons
    tablinks = document.getElementsByClassName("tablink");
    for (i = 0; i < tablinks.length; i++) {
      tablinks[i].classList.remove("selected");
    }

    // Show the specific tab content
    document.getElementById(pageName).style.display = "block";

    // Add the specific color to the button used to open the tab content
    elmnt.classList.add("selected");
  }

  /**
   * @description if the menu is shown, hide it. If it is invisible, show it!
   */
  static toggle() {
    //console.log('toggle');
    if (Menu.isShown()) {
      Menu.hide();
    }
    else {
      Menu.show();
    }
  }


  static hide() {
    menu.classList.remove("menuShow");
    menu.classList.add("menuHide");
  }


  static show() {
    menu.classList.add("menuShow");
    menu.classList.remove("menuHide");
  }

  static isShown() { return menu.classList.contains("menuShow"); }
}

/**
  * This class implements the switch between whiteboard and blackboard
  */
class BlackVSWhiteBoard {

  static init() {
    document.getElementById("canvasBackground").style.backgroundColor = "black";
    document.getElementById("whiteBoardSwitch").onclick = BlackVSWhiteBoard.switch;
    document.getElementById("blackBoardSwitch").onclick = BlackVSWhiteBoard.switch;
  }


  /**
   * switch between whiteboard and blackboard
   */
  static switch() {
    const previousBackgroundColor = document.getElementById("canvasBackground").style.backgroundColor;
    const backgroundColor = previousBackgroundColor == "white" ? "black" : "white";

    document.getElementById(backgroundColor + "BoardSwitch").hidden = true;
    document.getElementById(previousBackgroundColor + "BoardSwitch").hidden = false;


    console.log("previous background color was " + previousBackgroundColor);
    console.log("switch to " + backgroundColor + "board")
    palette.switchBlackAndWhite();
    document.getElementById("canvasBackground").style.backgroundColor = backgroundColor;

    BlackVSWhiteBoard._invertCanvas();
    UI.set_assistants_color();
    // should change magnets background color, but does not work correctly
    let magnets = document.getElementsByClassName("magnet");
    for (let i = 0; i < magnets.length; i++) {
      UI.set_color_and_background(magnets[i], palette.colors[0],palette.colors[0]);
    }
    //console.log(magnets);
  }


  /**
   * @dsecription invert the colors of the canvas (black becomes white, white becomes black, red becomes...)
   */
  static _invertCanvas() {
    const canvas = document.getElementById('canvas');
    const context = canvas.getContext('2d');

    const imageData = context.getImageData(0, 0, canvas.width, canvas.height);
    const data = imageData.data;

    for (let i = 0; i < data.length; i += 4) {
      // red
      data[i] = 255 - data[i];
      // green
      data[i + 1] = 255 - data[i + 1];
      // blue
      data[i + 2] = 255 - data[i + 2];
    }

    // overwrite original image
    context.putImageData(imageData, 0, 0);
  }

}

class Welcome {
  static init() {
    document.getElementById("welcome").onclick = Welcome.hide;
  }

  static isShown() {
    return !welcome.classList.contains("menuHide");
  }

  static hide() {
    //console.trace();
    welcome.classList.remove("menuShow");
    welcome.classList.add("menuHide");
  }

}

class TouchScreen {

  /**
   * 
   * @param {*} element 
   * @description makes that an element with events for the mouse, can be used on a smartphone/tablet
   */
  static addTouchEvents(element) {
    element.ontouchstart = TouchScreen._touchHandler;
    element.ontouchmove = TouchScreen._touchHandler;
    element.ontouchend = TouchScreen._touchHandler;
    element.ontouchcancel = TouchScreen._touchHandler;
  }
  /**
   * 
   * @param {*} event 
   * @description converts a touch event into a mouse event
   */
  static _touchHandler(event) {
    var touches = event.changedTouches,
        first = touches[0],
        type = "";
    if (UI.readonly && first.target.id=='canvas'){
      if (event.type=="touchstart"){
	UI.canvas_lastx=first.clientX; UI.canvas_lasty=first.clientY;
      }
      if (event.type=="touchmove" || event.type=="touchend"){
	//alert(first.target.id);
	var dx=first.clientX-UI.canvas_lastx;
	var dy=first.clientY-UI.canvas_lasty;
	// alert('mouse '+UI.canvas_lastx+','+first.clientX);
	if (dx*dx>dy*dy){
	  UI.canvas_lasty=first.clientY;
	  event.preventDefault();
	  if (dx*dx>2500){
	    if (dx<0) BoardManager.right(); else BoardManager.left();
	    UI.canvas_lastx=first.clientX;
	  }
	}
	else {
	  UI.canvas_lastx=first.clientX;
	}
      }
      return;
    }
    switch (event.type) {
    case "touchstart": type = "mousedown"; break;
    case "touchmove": type = "mousemove"; break;
    case "touchend": type = "mouseup"; break;
    default: return;
    }

    // initMouseEvent(type, canBubble, cancelable, view, clickCount, 
    //                screenX, screenY, clientX, clientY, ctrlKey, 
    //                altKey, shiftKey, metaKey, button, relatedTarget);

    var simulatedEvent = document.createEvent("MouseEvent");
    simulatedEvent.initMouseEvent(type, true, true, window, 1,
				  first.screenX, first.screenY,
				  first.clientX, first.clientY, false,
				  false, false, false, 0/*left*/, null);
    first.target.dispatchEvent(simulatedEvent);
    if (type=="mousemove")
      event.preventDefault();
  }
}

const SERVERADRESS = 'ws://tableaunoir.irisa.fr:8080';


class Share {
  static ws = undefined;
  static id = undefined;
  static canWriteValueByDefault = true;
  /**
   * 
   * @param {*} a function f 
   * @description tries to connect to the server, when the connection is made, it executes f
   */
  static tryConnect(f) {
    if (Share.ws != undefined)
      return;

    Share.ws = new WebSocket(SERVERADRESS);
    Share.ws.binaryType = "arraybuffer";

    Share.ws.onerror = () => { ErrorMessage.show("Impossible to connect to the server.") };

    Share.ws.onopen = f;
    Share.ws.onmessage = (msg) => {
      console.log("I received the message: ");
      Share._treatReceivedMessage(JSON.parse(msg.data));
    };

  }
  static init() {
    document.getElementById("shareButton").onclick = () => {
      if (!Share.isShared())
	Share.share();
    };

    document.getElementById("joinButton").onclick = () => {
      window.open(window.location, "_self")
    }

    document.getElementById("shareInEverybodyWritesMode").onclick = Share.everybodyWritesMode;
    document.getElementById("shareInTeacherMode").onclick = Share.teacherMode;

    if (window.location.origin.indexOf("github") < 0)
      document.getElementById('ShareGithub').hidden = true;

    if (Share.isSharedURL()) {
      let tryJoin = () => {
	try {
	  Share.id = Share.getIDInSharedURL();
	  if (Share.id != null) {
	    Share.join(Share.id);
	    document.getElementById("shareUrl").value = document.location;
	  }
	}
	catch (e) {
	  Share.ws = undefined;
	  ErrorMessage.show("Impossible to connect to the server", e);
	}

      }

      Share.tryConnect(tryJoin);
    }



  }


  /**
   * @returns true iff the board is shared with others
   */
  static isShared() {
    return Share.id != undefined;
  }


  /**
   * @description tries to connect the server to make a shared board
   */
  static share() {
    try {
      Share.tryConnect(() => Share.send({ type: "share" }));

      document.getElementById("shareInfo").hidden = false;
      document.getElementById("join").hidden = true;

    }
    catch (e) {
      Share.ws = undefined;
      ErrorMessage.show("Impossible to connect to the server", e);
    }

  }

  /**
   * @returns true if the userID of the current user is the minimum of all participants
   */
  static isSmallestUserID() {
    let minkey = "zzzzzzzzzzzzzzzz";
    for (let key in users) {
      if (key < minkey)
	minkey = key;
    }
    return (user.userID == minkey);
  }


  /**
   * 
   * @param {*} msg as an object
   * @description treats the msg received from the server
   */
  static _treatReceivedMessage(msg) {
    if (msg.type != "fullCanvas" && msg.type != "magnets" && msg.type != "execute")
      console.log("Server -> me: " + JSON.stringify(msg));
    else
      console.log("Server -> me: " + msg.type);
    switch (msg.type) {
    case "id": Share._setTableauID(msg.id); break;
    case "youruserid": // "your name is ..."
      Share._setMyUserID(msg.userid);

      document.getElementById("shareAndJoin").hidden = true;
      document.getElementById("shareInfo").hidden = false;
      document.getElementById("shareMode").hidden = false;

      break;
    case "user": //there is an existing user
      console.log("existing user: ", msg.userid)
      if (msg.userid == user.userID)
	throw "oops... an already existing user has the same name than me";

      users[msg.userid] = new User();
      Share.updateUsers();
      break;

    case "join": //a new user joins the group
      console.log("a new user is joining: ", msg.userid)
      // the leader is the user with the smallest ID

      users[msg.userid] = new User();
      Share.updateUsers();

      if (Share.isSmallestUserID()) {
	canvas.toBlob((blob) => Share.sendFullCanvas(blob, msg.userid));
	Share.sendFullCanvas(msg.userid);
	Share.sendMagnets(msg.userid);
	Share.execute("setUserCanWrite", [msg.userid, Share.canWriteValueByDefault]);
      }

      break;
    case "leave": users[msg.userid].destroy(); delete users[msg.userid]; Share.updateUsers(); break;
    case "fullCanvas": BoardManager.loadWithoutSave(msg.data); break;
    case "magnets":
      document.getElementById("magnets").innerHTML = msg.magnets;
      MagnetManager._installMagnetsNoMsg();
      break;
    case "execute": eval("ShareEvent." + msg.event)(...msg.params);
    }
  }

  /**
   * 
   * @param {*} msg as an object
   * @description send the message to server
   * 
   */
  static send(msg) {
    msg.id = Share.id;
    this.ws.send(JSON.stringify(msg));
  }

  /**
   * 
   * @param {*} blob 
   * @param {*} to 
   * @description send the blob of the canvas to the user to
   */
  static sendFullCanvas(blob, to) {
    Share.send({ type: "fullCanvas", data: canvas.toDataURL(), to: to }); // at some point send the blob directly
  }


  static sendMagnets(to) {
    if (Share.isShared())
      Share.send({ type: "magnets", magnets: document.getElementById("magnets").innerHTML, to: to }); // send the html code for all the magnets
  }


  /**
   * 
   * @param {*} event, an event name (string), that is a method of the class ShareEvent
   * @param {*} params an array of parameters
   * @description executes the event with the params, that is execute the method event of the class ShareEvent
   * with the params. Then send a message to server that this event should be executed for the other users as well
   */
  static execute(event, params) {
    if (!UI.readonly) UI.scroll_controls();
    function adapt(obj) {
      if (obj instanceof MouseEvent) {
	let props = [//'target', 'clientX', 'clientY', 'layerX', 'layerY', 
	  'pressure', 'offsetX', 'offsetY'];
	props.forEach(prop => {
	  Object.defineProperty(obj, prop, {
	    value: obj[prop],
	    enumerable: true,
	    configurable: true
	  });
	});
      }

      return obj;
    }
    eval("ShareEvent." + event)(...params);
    if (Share.isShared())
      Share.send({ type: "execute", event: event, params: params.map((param) => adapt(param)) });
  }


  static _setTableauID(id) {
    Share.id = id;
    let newUrl = document.location.href + "?id=" + id;
    history.pushState({}, null, newUrl);
    document.getElementById("shareUrl").value = newUrl;

    //document.getElementById("canvas").toBlob((blob) => Share.sendFullCanvas(blob));

  }


  static _setMyUserID(userid) {
    for (let key in users) {
      if (users[key] == user)
	delete users[key];
    }

    users[userid] = user;
    user.setUserID(userid);
    Share.updateUsers();
  }




  static updateUsers() {
    let i = 0;
    let keys = [];
    for (var key in users) {
      i++;
      if (key == user.userID)
	keys.push(key + "(me)");
      else
	keys.push(key);

    }
    document.getElementById("users").innerHTML = i + " users: " + keys.join("   ");
  }



  static isSharedURL() {
    let params = (new URL(document.location)).searchParams;
    return params.get('id') != null;
  }


  static getTableauNoirID() {
    if (Share.isSharedURL()) {
      return Share.getIDInSharedURL();
    }
    else
      return "local";
  }


  static getIDInSharedURL() {
    let params = (new URL(document.location)).searchParams;
    return params.get('id');
  }




  static join(id) {
    Share.send({ type: "join", id: id });
  }



  static teacherMode() {
    Share.setCanWriteForAllExceptMeAndByDefault(false);
  }


  static setCanWriteForAllExceptMeAndByDefault(bool) {
    for (let userid in users) {
      if (users[userid] != user)
	Share.execute("setUserCanWrite", [userid, bool]);
    }
    Share.canWriteValueByDefault = bool;
    Share.execute("setUserCanWrite", [user.userID, true]);
  }

  static everybodyWritesMode() {
    Share.setCanWriteForAllExceptMeAndByDefault(true);
  }
}


class ShareEvent {
  static mousedown(userId, evt) {
    users[userId].mousedown(evt);
  }

  static mousemove(userId, evt) {
    if (users[userId] == undefined)
      console.log("why is " + userId + " not declared?");
    users[userId].mousemove(evt);
  }

  static mouseup(userId, evt) {
    users[userId].mouseup(evt);
  }

  static setCurrentColor(userId, color) {
    users[userId].setCurrentColor(color);
  }

  static switchErase(userId) {
    users[userId].switchErase();
  }

  static switchChalk(userId) {
    users[userId].switchChalk();
  }

  static setUserCanWrite(userId, bool) {
    users[userId].setCanWrite(bool);
  }

  static magnetMove(idMagnet, x, y) {
    const el = document.getElementById(idMagnet);
    el.style.top = y + "px";
    if (x<0) x=0;
    el.style.left = x + "px";
  }


  static magnetsClear() {
    MagnetManager.clearMagnet();
  }
}

class ErrorMessage {

  static show(msg, msgplus) {
    console.log("error: " + msg);
    if (msgplus) console.log(msgplus);
    document.getElementById("error").hidden = false;
    document.getElementById("error").innerHTML = msg;

    const hide = () => { document.getElementById("error").hidden = true };
    document.getElementById("error").onclick = hide;
    setInterval(hide, 5000);
  }
}
