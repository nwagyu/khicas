function $id(id) { return document.getElementById(id); }

var UI = {
  readonly:false, // set to true will disable tableau changes
  locked:false, // lock right and left
  scroll_controls:function(){
    //console.log('scroll');
    $id('controls').scrollIntoView();
  },
  set_warnlink:function(val){
    UI.warn_link=val;
    UI.createCookie('xcas_warn_link',val?'1':'0')
    $id('warnlink_yes').style.fontWeight=val?"900":"normal";
    $id('warnlink_no').style.fontWeight=val?"normal":"900";    
  },    
  set_readonly:function(val){
    UI.readonly=val;
    $id('addmagnetbutton').style.display=val?'none':'inline';
    $id('pc_controls_readonly_hidden').style.display=val?'none':'inline';
    $id('readonly_yes').style.fontWeight=val?"900":"normal";
    $id('readonly_no').style.fontWeight=val?"normal":"900";
    controls.hidden=false;
    UI.scroll_controls();
    if (val){
      $id("readonly_msg").style.display='inline';
      setTimeout( ()=>{$id("readonly_msg").style.display='none';},10000);
    }
    else
      $id('readonly_msg').style.display='none';      
  },
  gr2d_ncanvas:0,
  disable3d:0, // -1==3d not supported
  tableaucm:0, // 1 if codemirror focused
  debug:0,
  Datestart:Date.now(),
  ready: false,
  colorinit:false,
  focusaftereval: false,
  docprefix: "https://www-fourier.univ-grenoble-alpes.fr/%7eparisse/giac/doc/fr/cascmd_fr/",
  mailto: '',
  from: '',
  base_url: "https://www-fourier.univ-grenoble-alpes.fr/%7eparisse/",
  //forum_url: "http://xcas.e.univ-grenoble-alpes.fr/XCAS/viewforum.php?f=25",
  forum_url: "https://xcas.univ-grenoble-alpes.fr/forum/viewforum.php?f=25",
  // forum_url: "http://xcas.e.univ-grenoble-alpes.fr/XCAS/posting.php?mode=post&f=12&subject=session",
  //forum_url: "http://xcas.univ-grenoble-alpes.fr/forum/posting.php?mode=post&f=12&subject=session",
  forum_warn: true,
  hints: false, // set this to true/false in show-hint.js
  createCookie: function (name, value, days) {
    if (window.localStorage) return localStorage.setItem(name, value);
    if (days) {
      var date = new Date();
      date.setTime(date.getTime() + (days * 24 * 60 * 60 * 1000));
      var expires = "; expires=" + date.toGMTString();
    }
    else var expires = "";
    document.cookie = name + "=" + value + expires + "; path=/";
  },
  readCookie: function (name) {
    if (window.localStorage) {
      var tmp = localStorage.getItem(name);
      if (tmp != null) return tmp;
    }
    var nameEQ = name + "=";
    var ca = document.cookie.split(';');
    for (var i = 0; i < ca.length; i++) {
      var c = ca[i];
      while (c.charAt(0) == ' ') c = c.substring(1, c.length);
      if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length, c.length);
    }
    return null;
  },
  eraseCookie: function (name) {
    if (window.localStorage) return localStorage.removeItem(name);
    createCookie(name, "", -1);
  },
  compactlink:true,
  palette_colors:["white", "yellow", "orange", "rgb(100, 172, 255)", "Crimson", "Plum", "LimeGreen"],
  point_in_polygon:function(testx,testy,vertx,verty){
    // is testx,testy inside polygon of vertices coordinates in lists vertx verty
    // (https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html)
    let i, j, c = 0,nvert=vertx.length;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
      if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	   (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
	c = !c;
    }
    return c;
  },
  cut_stack:function(xmin,xmax,ymin,ymax,points=0){
    // remove delineations starting in points or in xmin..xmax, ymin..ymax
    // returns the delineations
    let tab=BoardManager.cancelStack.stack;
    //console.log('cut_stack',xmin,xmax,ymin,ymax,tab);
    let delineations='';
    if (xmin>xmax){ [xmin,xmax]=[xmax,xmin];}
    if (ymin>ymax){ [ymin,ymax]=[ymax,ymin];}
    for (let i=0;i<tab.length;++i){
      if (!tab[i])
	continue;
      let cur=UI.delineation_start(tab[i]);
      if (points?
	  Delineation.inPolygon({x:cur[0],y:cur[1]},points):
	  (cur[0]>=xmin && cur[0]<xmax && cur[1]>=ymin && cur[1]<=ymax)
	 ){
	// remove from stack, put in delineations
	delineations += tab[i]+'&';
	tab.splice(i,1); --i; if (BoardManager.cancelStack.currentIndex>i) BoardManager.cancelStack.currentIndex--; if (BoardManager.cancelStack.n>i) BoardManager.cancelStack.n--;
      }
    }
    return delineations;
  },
  console_show_stack:function(){
    console.log(BoardManager.cancelStack.stack.length,BoardManager.cancelStack.n,BoardManager.cancelStack.currentIndex,BoardManager.cancelStack.stack);
  },
  delineation_start:function(s,DX=0,DY=0){ // return [startX,startY]
    let l=s.length,pos=s.search(','),lastX=0,lastY=0;
    if (s[0]>'9'){ // skip color
      s=s.substr(pos+1,l-pos-1);
      pos=s.search(',');
    }
    lastX=UI.parse_int(s.substr(0,pos))+DX;
    s=s.substr(pos+1,l-pos-1);
    pos=s.search(',');
    lastY=UI.parse_int(s.substr(0,pos))+DY;
    return [lastX,lastY];
  },
  decode_polygon:function(s,DX=0,DY=0,truncate=0){
    // s is a string of format [color_if_not_white,]startX,startY,
    // followed by points of the delineation in a "compact" format [pressure]dxdy
    // where dx and dy are encoded with UI.compact
    // returns lx,ly,lp,color,s where
    // lx, ly is the list of x,y stroke coordinate
    // lp list of pressures
    // color is the color
    // s is the part of initial string after color+initial absolute coordinates
    let lx=[],ly=[],lp=[];
    let l=s.length,pos=s.search(','),lastX=0,lastY=0;
    let color="white",pressure=0.5;
    if (s[0]>'9'){
      color=s.substr(0,pos);
      s=s.substr(pos+1,l-pos-1);
      pos=s.search(',');
    }
    lastX=UI.parse_int(s.substr(0,pos))+DX;
    s=s.substr(pos+1,l-pos-1);
    pos=s.search(',');
    lastY=UI.parse_int(s.substr(0,pos))+DY;
    lx.push(lastX); ly.push(lastY); lp.push(pressure);
    // console.log(color,lastX,lastY); return;
    l=s.length; s=s.substr(pos+1,l-pos-1);
    // now s has format [digit]majuscules minuscules
    l=s.length;
    for (pos=0;pos<l;){
      let n=s.charCodeAt(pos);
      if (n>=48 && n<=57){
	pressure=(n-48)/10+0.05;
	//console.log('pressure',pressure);
	++pos;
      }
      let dx=0,dy=0;
      for (;pos<l;++pos){
	n=s.charCodeAt(pos);
	if (n>96)
	  break;
	if (n<65){
	  console.log('invalid char');
	  return [lx,ly,lp,color];
	}
	dx *= 13;
	dx += (n-77); // 65+12 where 'A'==65
      }
      // console.log(dx);
      for (;pos<l;++pos){
	n=s.charCodeAt(pos);
	if (n<97)
	  break;
	if (n>122){
	  console.log('invalid char');
	  return [lx,ly,lp,color];
	}
	dy *= 13;
	dy += (n-109); // 97+12 where 'A'==65
      }
      // console.log(dy);
      let curX=lastX+dx;
      let curY=lastY+dy;
      //console.log(dx,dy);
      lastX=curX; lastY=curY;
      lx.push(lastX); ly.push(lastY); lp.push(pressure);
      if (truncate && lx.length>=truncate)
	break;
    }
    return [lx,ly,lp,color,s];
  },
  find_whxy:function(tab){
    // find canvas size that contains these relative handwritings
    // and dx,dy shift to write s from 0,0
    let xmin=1e307,xmax=-1e307,ymin=1e307,ymax=-1e307,lastX,lastY;
    for (let i=0;i<tab.length;++i){
      let s=tab[i];
      // console.log('handwrite_rectangle',i,s);
      let l=s.length,pos=s.search(','),lastX=0,lastY=0;
      if (s[0]>'9'){
	s=s.substr(pos+1,l-pos-1);
	pos=s.search(',');
      }
      lastX=UI.parse_int(s.substr(0,pos));
      s=s.substr(pos+1,l-pos-1);
      pos=s.search(',');
      lastY=UI.parse_int(s.substr(0,pos));
      if (xmin>lastX)
	xmin=lastX;
      if (xmax<lastX)
	xmax=lastX;
      if (ymin>lastY)
	ymin=lastY;
      if (ymax<lastY)
	ymax=lastY;
      // console.log(lastX,lastY); // return;
      l=s.length; s=s.substr(pos+1,l-pos-1);
      // now s has format [digit]majuscules minuscules
      l=s.length;
      for (pos=0;pos<l;){
	let n=s.charCodeAt(pos);
	if (n>=48 && n<=57){
	  //console.log('pressure',pressure);
	  ++pos;
	}
	let dx=0,dy=0;
	for (;pos<l;++pos){
	  n=s.charCodeAt(pos);
	  if (n>96)
	    break;
	  if (n<65){
	    console.log('invalid char');
	    return c;
	  }
	  dx *= 13;
	  dx += (n-77); // 65+12 where 'A'==65
	}
	// console.log(dx);
	for (;pos<l;++pos){
	  n=s.charCodeAt(pos);
	  if (n<97)
	    break;
	  if (n>122){
	    console.log('invalid char');
	    return c;
	  }
	  dy *= 13;
	  dy += (n-109); // 97+12 where 'A'==65
	}
	// console.log(dy);
	let curX=lastX+dx;
	let curY=lastY+dy;
	//console.log(dx,dy);
	lastX=curX; lastY=curY;
	if (xmin>lastX)
	  xmin=lastX;
	if (xmax<lastX)
	  xmax=lastX;
	if (ymin>lastY)
	  ymin=lastY;
	if (ymax<lastY)
	  ymax=lastY;
      }
    }
    // console.log(xmin,xmax,ymin,ymax);
    w=xmax+1-xmin;
    h=ymax+1-ymin;
    dx=-xmin;
    dy=-ymin;
    return [w,h,dx,dy];
  },
  handwrite:function(s_in,w=0,h=0,dx=0,dy=0){
    // console.log('handwrite',s_in);
    // return a canvas from handwritings encoded in s
    let c=document.createElement("CANVAS");
    c.onfocus=(e) => { this.blur(); this.parentNode.focus();};
    if (s_in.length && s_in[s_in.length-1]=='&')
      s_in=s_in.substr(0,s_in.length-1);
    let tab=s_in.split('&');
    /*
    for (let i=0;i<tab.length;++i){
      let s_=UI.parse_util(tab[i]);
      console.log(tab[i],s_);
      tab[i]=s_.substr(1,s_.length-1); // remove '=' at start
    }
    */
    if (w==0 || h==0)
      [w,h,dx,dy]=UI.find_whxy(tab);
    c.width=w; c.height=h;
    let ctx=c.getContext("2d");
    for (let i=0;i<tab.length;++i)
      UI.draw_encoded(tab[i],ctx,dx,dy);
    return c;
  },
  handwritings_string:function(){ // return handwriting encoded in a string 
    let s='';
    let n=BoardManager.cancelStack.currentIndex;
    let data=BoardManager.cancelStack.stack;
    for (let i=0;i<n;++i){
      if (UI.compactlink){
	// compact format for color, position
	let cs=data[i];
	if (cs){
	  //console.log(cs);
	  let l=cs.length,pos=cs.search(','),lastX=0,lastY=0;
	  let colorn='0',pressure=0.5;
	  if (cs[0]>'9'){
	    let color=cs.substr(0,pos);
	    for (let j=0;j<UI.palette_colors.length;++j){
	      if (color==UI.palette_colors[j]){
		colorn=j;
		break;
	      }
	    }
	    cs=cs.substr(pos+1,l-pos-1);
	    pos=cs.search(',');
	  }
	  lastX=UI.parse_int(cs.substr(0,pos));
	  cs=cs.substr(pos+1,l-pos-1);
	  pos=cs.search(',');
	  lastY=UI.parse_int(cs.substr(0,pos));
	  // console.log(color,lastX,lastY); return;
	  l=cs.length; cs=cs.substr(pos+1,l-pos-1);
	  s += colorn+UI.compact(lastX)+UI.compact(lastY)+cs+'&';
	}
      }
      else
	s += '='+data[i]+'&';
    }
    return s;
  },
  warn_link:true,
  magnetlink:function(magnets){
    let s='';
    for (let i = 0; i < magnets.length; i++) {
      let cur=magnets[i];
      if (typeof cur=='string'){
	s += 'handwriting=0,0,'+encodeURIComponent(cur)+'&';
	continue;
      }
      if (cur.style.visibility=='hidden')
	continue;
      let bg='';
      if ( cur.style.backgroundColor!='' && cur.style.backgroundColor!=((palette.colors[0]== "white") ? "rgba(0, 0, 0, 0.9)" : "rgba(1, 1, 1, 0.9)") ){
	bg=cur.style.backgroundColor+',';
      }
      let mx=parseInt(cur.style.left)+container.scrollLeft,my=parseInt(cur.style.top);
      if (cur.tagName=='IMG'){
	//console.log(cur);
	s += 'img='+bg+mx+','+my+','+encodeURIComponent(cur.src)+'&';
	continue;
      }
      if (cur.tagName=='DIV' && cur.hasChildNodes()){
	if (cur.firstChild.tagName=='IMG'){
	  //console.log(cur);
	  s += 'img='+bg+mx+','+my+','+encodeURIComponent(cur.firstChild.src)+'&';
	  continue;
	}
	if (cur.firstChild.tagName=='svg'){
	  s += 'svg='+bg+mx+','+my+','+encodeURIComponent(cur.innerHTML)+'&';
	  continue;
	}
	if (cur.firstChild.nextSibling && cur.firstChild.nextSibling.tagName=='CANVAS'){
	  let adds='handwriting='+bg+mx+','+my+','+encodeURIComponent(cur.firstChild.innerText)+'&';
	  s += adds;
	  //console.log('link hand i bg mx my',i,bg,mx,my,adds);
	  continue;
	}
	if (cur.firstChild.nextSibling && cur.firstChild.nextSibling.type=='textarea'){
	  let t=cur.firstChild.innerHTML;
	  let v=encodeURIComponent(cur.firstChild.nextSibling.value);
	  v = bg+mx+','+my+','+v;
	  if (t=='&gt;')
	    v = 'cas='+v;
	  if (t=='&gt;&gt;')
	    v = 'js='+v;
	  if (t=='&gt;&gt;&gt;')
	    v = 'py='+v;
	  s += v+'&';
	  // console.log('save eval magnet:',i,v);
	  continue;
	}
	if (cur.firstChild.type=='textarea'){
	  let v=encodeURIComponent(cur.firstChild.value);
	  v = bg+mx+','+my+','+v;
	  s += 'textarea='+v+'&';
	  continue;
	}
	cur=cur.firstChild;
	if (cur.tagName=='SPAN' && cur.hasChildNodes()){
	  //console.log(cur.firstChild);
	  if (cur.firstChild.type=='textarea'){
	    let v=encodeURIComponent(cur.firstChild.value);
	    //console.log(cur.parentNode.style.color,cur);
	    if (bg!='' && cur.parentNode.style.color!='')
	      bg += cur.parentNode.style.color+',';
	    v = bg+mx+','+my+','+v;
	    s += 'comment='+v+'&';
	    continue;
	  }
	}
      }
    }
    return s;
  },
  magnetize_all:function(){
    // create an ordered list of clones of visible magnets and handwritings delineations
    let magnets = document.getElementsByClassName("magnet");
    let tab=[],tabhw=[''];
    let largeur=BoardManager.scrollQuantity();
    for (let i = 0; i < magnets.length; i++) {
      if (magnets[i].style.visibility=='hidden')
	continue;
      tab.push(magnets[i]);
      tabhw.push(''); // empty handwritings after i
    }
    let N=tab.length;
    tab.sort(UI.order_magnets);
    // console.log('magnetize_all',tab);
    // now look in cancelstack and add handw in tabhw at the right position
    let hw=BoardManager.cancelStack.stack;
    for (let i=0;i<hw.length;++i){
      let s=hw[i];
      if (N==0){
	tabhw[0] += s+'&';
	continue;
      }
      let lxypc=UI.decode_polygon(s,0,0,2);
      let lx=lxypc[0],ly=lxypc[1];
      if (lx.length<1 || ly.length!=lx.length)
	continue;
      if (lx.length==2 && lx[0]==lx[1] && ly[1]-ly[0]>largeur) // vertical lines
	continue;
      lx=lx[0]; ly=ly[0];
      //console.log('magnetize_all',s,lx,ly);
      lx=Math.floor(lx/largeur)*largeur;
      let j=0;
      for (;j<tab.length;++j){ // dichotomy would be faster
	let lX=parseInt(tab[j].style.left)+parseInt(container.scrollLeft);
	lX=Math.floor(lX/largeur)*largeur;
	// console.log(lx,ly,j,lX,parseInt(tab[j].style.top));
	if (lx<lX || (lx==lX && ly<parseInt(tab[j].style.top))){
	  tabhw[j] += s+'&';
	  break;
	}
      }
      if (j==tab.length)
	tabhw[j] += s+'&';
    }
    // merge with tab
    let mergetab=[]; 
    for (let i=0;i<tab.length;++i){
      if (tabhw[i]!='')
	mergetab.push(tabhw[i]);
      mergetab.push(tab[i]);
    }
    if (tabhw[tab.length]!='')
      mergetab.push(tabhw[tab.length]);
    return mergetab;
  },
  makelink: function (convert) { // convert=1 if all handwritings are converted
    if (UI.warn_link && !UI.readonly){
      UI.warn_link=false;
      if (!confirm(UI.langue==-1?"Lien de partage mail/forum créé. Attention l'envoi par mail peut faire croire à du spam. Conserver cet avertissement ?":"Mail/forum share link created. Beware that sending by mail might be considered to be spam. Keep this warning?"))
	UI.createCookie('xcas_warn_link','0');
    }
    //console.log('makelink',start);
    $id("share_submenu").style.display='inline';
    let s = 'python=';
    if (UI.python_mode) s += (UI.python_mode+'&'); else s += (UI.micropy+'&');;
    s += 'radian='+(UI.radian_mode?1:0)+'&';
    s += 'width='+canvas.width+'&';
    s += 'height='+canvas.height+'&';
    s += 'blackboard='+palette.colors[0]+'&';
    if (convert){
      s += UI.magnetlink(UI.magnetize_all());
    }
    else {
      s += UI.handwritings_string();
      // add magnets
      s += UI.magnetlink($id("magnets").children); //document.getElementsByClassName("magnet");
    }
    //console.log(s);
    return s;
  },
  link: function (start) {
    let s = UI.makelink(0),sxcas=UI.makelink(1);
    //if ($id('variables').style.display != 'none') UI.listvars(3);
    //console.log('link',sxcas);
    UI.createCookie('xcas_tableau', s, 365);
    if (s.length > 0) {
      let s2 = "#exec&" + s;
      let smail;
      if (UI.langue == -1)
        smail = UI.base_url + "tableaufr.html#exec&";
      else
        smail = UI.base_url + "tableauen.html#exec&";
      smail = smail + s;
      let sforum;
      if (UI.langue == -1) {
        sforum = UI.base_url + "tableaufr.html#exec&" + s;
	sxcas = "xcasfr.html#" + sxcas;
        s = UI.base_url + "tableaufr.html#" + s;
      }
      else {
        sforum = UI.base_url + "tableauen.html#exec&" + s;
	sxcas = "xcasen.html#" + sxcas;
        s = UI.base_url + "tableauen.html#" + s;
      }
      if (0) sxcas = UI.base_url + sxcas;
      sxcas = '<a href="'+ sxcas+'" target="_blank">Xcas</a>';
      //s=encodeURIComponent(s); // does not work innerHTML will add a prefix
      //var sforum=encodeURIComponent('[url]'+s+'[/url]');
      sforum = '[url=' + sforum + ']tableau[/url]';
      //console.log(sforum);
      $id('theforumlink').innerHTML = sforum;
      var copy = "<button class='spanlink' title=";
      copy += UI.langue == -1 ? "'Partager cette session sur le forum'" : "'Share this session on the forum'";
      copy += " onclick='var tmp=$id(\"theforumlink\"); tmp.style.display=\"inline\";tmp.select();document.execCommand(\"copy\");tmp.style.display=\"none\"; ";
      if (UI.forum_warn) {
        UI.forum_warn = false;
        copy += UI.langue == -1 ? "alert(\"Le lien de la session a été copié dans le presse-papier\");" : "alert(\"Clipboard contains a link to session\");";
      }
      copy += "var win=window.open(\"" + UI.forum_url + "\", \"_blank\");'>Forum</button>";
      //console.log(copy);
      if (window.location.href.substr(0, 4) == 'file' && !UI.detectmob()) {
        $id('thelink').innerHTML = '<a href="#" title="Nouveau tableau" target="_blank">+</a> <a title="Clone session" href="' + s + '" target="_blank">x2</a> <a title="Local clone" href="' + s2 + '" target="_blank">local</a>' + copy+sxcas;//+',<a href="http://xcas.e.univ-grenoble-alpes.fr/XCAS/posting.php?mode=post&f=12&subject=session&message='+encodeURIComponent(sforum)+'" target="_blank">forum</a>,';
      }
      else
        $id('thelink').innerHTML = '<a href="' + s + '" target="_blank">x2</a> ' + copy+sxcas;
      var mailurl;
      if (UI.from.length > 9 && UI.from.substr(UI.from.length - 9, 9) == "gmail.com")
        mailurl = 'https://mail.google.com/mail/?view=cm&fs=1&tf=1&source=mailto&su=session+Xcas&to=' + UI.mailto;
      else
        mailurl = 'mailto:' + UI.mailto + '?subject=tableaunoirXcas';
      mailurl += '&body=Bonjour%0d%0aVeuillez suivre ce lien : <' + encodeURIComponent(smail) + '>';
      mailurl = '<a title="E-mail session" href="' + mailurl + '" target="_blank"> &#x2709; </a>';
      $id('themailto').innerHTML = mailurl;
      $id('themailtopc').innerHTML = mailurl;
      $id('themailtomenu').innerHTML = mailurl;
      $id('thelinkmenu').innerHTML = $id('thelink').innerHTML;
      $id('themailtopc').style.display='inline';
      $id('thelink').style.display='inline';
    }
  },
  toggle:function(f){
    if (f.style.display=='none')
      f.style.display='block';
    else {
      f.style.display='none';
      if (UI.focused!=0){ UI.focused.focus(); UI.scroll_controls();}
    }
  },
  menu_modules:function(){
    UI.close_menus();$id('assistant_modules').style.display='block';
  },
  assistant_list:['online_help','assistant_cas','assistant_menu','assistant_math','assistant_arit','assistant_proba','assistant_rewrite','assistant_calculus','assistant_graph','assistant_prog','assistant_plotfunc','assistant_plotfunc1var','assistant_plotfunc2var','assistant_plotparam','assistant_plotparam1var','assistant_plotparam2var','assistant_plotimplicit','assistant_plotfield','assistant_statgraph','assistant_modules','assistant_math_module','assistant_arit_module','assistant_linalg_module','assistant_turtle_module','assistant_graphic_module','assistant_series','assistant_limit','assistant_int','assistant_diff','assistant_sum','assistant_solve','assistant_rsolve','assistant_rsolve','assistant_suites','assistant_fixe','assistant_desolve','assistant_test','assistant_boucle','assistant_pour','assistant_tantque','assistant_def','assistant_seq','assistant_tabvar','assistant_tabvarfunc','assistant_tabvarparam'],
  set_field_color:function(fieldid,color){
    let f=$id(fieldid);
    if (f){
      f.style.color=color;
    }
  },
  change_magnet_color:function(magnet){
    let color=palette.colors[0];
    UI.set_color_and_background(magnet, color,nextBackgroundColor(magnet.style.backgroundColor));
  },    
  set_color_and_background:function(field,color,back){
    //console.log(color,back);
    if (back==color || (color=="white" && back=='rgba(1, 1, 1, 0.1)') ||
	(color=="black" && back=='rgba(0, 0, 0, 0.9)') ){
      if (color=="white") color='black';
      else if (color=="black") color='white';
    }
    //console.log('after',color,back);
    field.style.color=color;
    field.style.backgroundColor=back;
    if (field.hasChildNodes()){
      for (let o of field.children)
	UI.set_color_and_background(o,color,back);
    }
  },
  set_assistants_color:function(){
    let color=$id('canvas').backgroundColor=='white'?'black':'white';
    console.log('set_assistant_color',color);
    for (let j=0;j<UI.assistant_list.length; ++j)
      UI.set_field_color(UI.assistant_list[j],color);
    UI.link();
  },
  hide_field:function(fieldid){
    let f=$id(fieldid);
    if (f) f.style.display='none';
  },
  close_menus:function(){
    //console.trace();
    for (let j=0;j<UI.assistant_list.length; ++j)
      UI.hide_field(UI.assistant_list[j]);
    //setTimeout(UI.scroll_controls(),100);
  },
  open_menu:function(toggle=true){
    $id('assistant_keys').style.display='block';
    if (UI.micropy){
      if (toggle && $id('assistant_prog').style.display!='none')
	$id('assistant_prog').style.display='none';
      else
	$id('assistant_prog').style.display='block';
    }
    else {
      if (toggle && $id('assistant_cas').style.display!='none'){
	$id('assistant_cas').style.display='none';
      }
      else {
	$id('assistant_cas').style.display='block';
	$id('assistant_menu').style.display='block';
      }
    }
  },
  focused: 0,
  savefocused: 0,
  restore_focus: function(){
    console.log('restore_focus',UI.savefocused);
    UI.focused=UI.savefocused;
    UI.savefocused=0;
    try {
      UI.focused.focus();
    }
    catch (err){}
  },
  clear_focus:function(){
    UI.focused=0;
    UI.close_menus();
    // $id('assistant_cas').style.display='none';
    $id('assistant_keys').style.display='none';
  },
  selectionne: function () {
    if (!UI.focused) return;
    UI.focused.focus();
    if (UI.focused.type != 'textarea') {
      if (UI.focused.execCommand)
        UI.focused.execCommand('selectAll');
    }
    else {
      UI.focused.select();
      UI.selection = UI.focused.value;
    }
  },
  set_focus: function (s) {
    console.log('set_focus',s);
    UI.savefocused=UI.focused;
    var tmp = $id(s);
    UI.focused = tmp;
    UI.selectionne(); // if (!UI.is_touch_device()) tmp.focus();
  },
  assistant_pour_ok: function () {
    UI.restore_focus();
    var st = $id('pourvarstep').value;
    if (UI.python_mode) {
      var sup=eval($id('pourvarmax').value)+1;
      UI.insert(UI.focused, '\nfor ' + $id('pourvarname').value + ' in range(' + $id('pourvarmin').value + ',' + sup);
      if (st.length) st = ',' + st;
      UI.insert(UI.focused, st + '):');
      UI.indentline(UI.focused);
      UI.insert(UI.focused, '\n');
      UI.indentline(UI.focused);
    }
    else {
      var tmp = '\npour ' + $id('pourvarname').value + ' de ' + $id('pourvarmin').value + ' jusque ' + $id('pourvarmax').value;
      //console.log(tmp);
      UI.insert(UI.focused, tmp);
      UI.indentline(UI.focused);
      if (st.length) UI.insert(UI.focused, ' pas ' + st);
      UI.insert(UI.focused, ' faire\n\nfpour;');
      UI.indentline(UI.focused);
      UI.moveCaretUpDown(UI.focused, -1);
      UI.indentline(UI.focused);
    }
    UI.focused.focus();
    $id('assistant_pour').style.display = 'none';
    $id('assistant_boucle').style.display = 'none';
  },
  assistant_tantque_ok: function () {
    UI.restore_focus();
    if (UI.python_mode) {
      UI.insert(UI.focused, '\nwhile ' + $id('tantquecond').value + ':');
      UI.indentline(UI.focused);
      UI.insert(UI.focused, '\n');
    }
    else {
      UI.insert(UI.focused, '\ntantque ' + $id('tantquecond').value + ' faire\n\nftantque;');
      UI.indentline(UI.focused);
      UI.moveCaretUpDown(UI.focused, -1);
      UI.indentline(UI.focused);
      UI.moveCaretUpDown(UI.focused, -1);
      UI.indentline(UI.focused);
      UI.moveCaretUpDown(UI.focused, 1);
    }
    UI.focused.focus();
    $id('assistant_tantque').style.display = 'none';
    $id('assistant_boucle').style.display = 'none';
  },
  assistant_test_ok: function () {
    UI.restore_focus();
    if (UI.python_mode) {
      UI.insert(UI.focused, '\nif ' + $id('sicond').value + ':');
      UI.indentline(UI.focused);
      UI.insert(UI.focused, '\n' + $id('sialors').value);
    } else UI.insert(UI.focused, '\nsi ' + $id('sicond').value + ' alors ' + $id('sialors').value);
    UI.indentline(UI.focused);
    var tmp = $id('sisinon').value;
    if (tmp.length) {
      if (UI.python_mode) {
        UI.indentline(UI.focused);
        UI.insert(UI.focused, '\nelse:')
        UI.indentline(UI.focused); // should remove 2 spaces at start
        UI.insert(UI.focused, '\n' + tmp);
      } else {
        UI.insert(UI.focused, ' sinon ' + tmp);
        UI.insert(UI.focused, ' fsi;\n');
      }
    }
    else {
      if (!UI.python_mode)
        UI.insert(UI.focused, ' fsi;\n');
    }
    UI.indentline(UI.focused);
    UI.focused.focus();
    $id('assistant_test').style.display = 'none';
  },
  assistant_def_ok: function () {
    UI.restore_focus();
    var loc = $id('localvars').value;
    var fc = $id('funcname').value;
    var argu = $id('argsname').value;
    var ret = $id('returnedvar').value;
    if (UI.python_mode) {
      UI.insert(UI.focused, 'def ' + fc + '(' + argu + '):');
      UI.indentline(UI.focused);
      if (0 && loc.length != 0) {
        UI.insert(UI.focused, '\n# local ' + loc);
        UI.indentline(UI.focused);
        UI.insert(UI.focused, '\n\nreturn ' + ret);
        UI.indentline(UI.focused);
        UI.moveCaretUpDown(UI.focused, -1);
        UI.indentline(UI.focused);
      }
      else {
        UI.insert(UI.focused, '\nreturn ' + ret);
        UI.indentline(UI.focused);
      }
    }
    else {
      if (loc.length == 0)
        UI.insert(UI.focused, 'fonction ' + fc + '(' + argu + ')\n \nffonction:;\n'); // was fc + '(' + argu + '):=' + ret + ';');
      else {
        UI.insert(UI.focused, 'fonction ' + fc + '(' + argu + ')\n  local ' + loc + ';\n  \n  retourne ' + ret + ';\nffonction:;\n');
        UI.moveCaretUpDown(UI.focused, -3);
        UI.moveCaret(UI.focused, 2);
      }
    }
    $id('assistant_def').style.display = 'none';
    UI.focused.focus();
  },
  assistant_seq_ok: function () {
    UI.restore_focus();
    var tmp = 'seq(' + $id('seqexpr').value + ',' + $id('seqvarname').value + ',' + $id('seqvarmin').value + ',' + $id('seqvarmax').value;
    UI.insert(UI.focused, tmp);
    tmp = $id('seqvarstep').value;
    if (tmp.length) UI.insert(UI.focused, ',' + tmp);
    UI.insert(UI.focused, ')');
    $id('assistant_seq').style.display = 'none';
  },
  assistant_series_ok: function () {
    UI.restore_focus();
    var tmp = 'series(' + $id('seriesexpr').value + ',' + $id('seriesvarname').value + '===' + $id('seriesvarlim').value + ',' + $id('seriesvarorder').value;
    UI.insert(UI.focused, tmp);
    tmp = $id('seriesvarstep').value;
    if (tmp.length) UI.insert(UI.focused, ',' + tmp);
    UI.insert(UI.focused, ')');
    $id('assistant_series').style.display = 'none';
  },
  assistant_limit_ok: function () {
    UI.restore_focus();
    var tmp = 'limit(' + $id('limitexpr').value + ',' + $id('limitvarname').value + ',' + $id('limitvarlim').value;
    UI.insert(UI.focused, tmp);
    tmp = $id('limitvardir').value;
    if (tmp.length) UI.insert(UI.focused, ',' + tmp);
    UI.insert(UI.focused, ')');
    $id('assistant_limit').style.display = 'none';
  },
  assistant_int_ok: function () {
    UI.restore_focus();
    var tmp = $id('intexpr').value;
    if (tmp.length) {
      tmp = 'integrate(' + tmp + ',' + $id('intvarname').value;
      UI.insert(UI.focused, tmp);
      tmp = $id('intvarmin').value;
      if (tmp.length) UI.insert(UI.focused, ',' + tmp + ',' + $id('intvarmax').value);
      UI.insert(UI.focused, ')');
    } else UI.insert(UI.focused, 'integrate(');
    $id('assistant_int').style.display = 'none';
  },
  assistant_diff_ok: function () {
    UI.restore_focus();
    var tmp = $id('diffexpr').value;
    if (tmp.length) {
      tmp = 'diff(' + tmp + ',' + $id('diffvarname').value;
      UI.insert(UI.focused, tmp);
      tmp = $id('diffnumber').value;
      if (tmp.length) UI.insert(UI.focused, ',' + tmp);
      UI.insert(UI.focused, ')');
    } else UI.insert(UI.focused, 'diff(');
    $id('assistant_diff').style.display = 'none';
  },
  assistant_sum_ok: function () {
    UI.restore_focus();
    var tmp = $id('sumexpr').value;
    if (tmp.length) {
      tmp = 'sum(' + tmp + ',' + $id('sumvarname').value;
      UI.insert(UI.focused, tmp);
      tmp = $id('sumvarmin').value;
      if (tmp.length) UI.insert(UI.focused, ',' + tmp + ',' + $id('sumvarmax').value);
      UI.insert(UI.focused, ')');
    } else UI.insert(UI.focused, 'sum(');
    $id('assistant_sum').style.display = 'none';
  },
  assistant_solve_ok: function () {
    UI.restore_focus();
    if ($id('solveC').style.display == 'inline') UI.insert(UI.focused, 'c');
    if ($id('solvenum').style.display == 'inline') UI.insert(UI.focused, 'f');
    UI.insert(UI.focused, 'solve(' + $id('solveeq').value);
    UI.indentline(UI.focused);
    var tmp = $id('solvevar').value;
    if (tmp.length) UI.insert(UI.focused, ',' + tmp);
    UI.insert(UI.focused, ')');
    $id('assistant_solve').style.display = 'none';
  },
  assistant_rsolve_ok: function () {
    UI.restore_focus();
    UI.insert(UI.focused, 'rsolve(' + $id('rsolveeq').value + ',' + $id('rsolvevar').value);
    var tmp = $id('rsolveinit').value;
    if (tmp.length) UI.insert(UI.focused, ',[' + tmp + ']');
    UI.insert(UI.focused, ')');
    $id('assistant_rsolve').style.display = 'none';
    $id('assistant_suites').style.display = 'none';
  },
  assistant_desolve_ok: function () {
    UI.restore_focus();
    var tmpeq = $id('desolveeq').value;
    var tmpt = $id('desolvevar').value;
    var tmpy = $id('desolvey').value;
    var tmp = $id('desolveinit').value;
    if (tmp.length) tmp = 'desolve([' + tmpeq + ',' + tmp + ']'; else tmp = 'desolve(' + tmpeq;
    tmp += ',' + tmpt + ',' + tmpy + ')';
    UI.insert(UI.focused, tmp);
    $id('assistant_desolve').style.display = 'none';
  },
  assistant_plotfunc1var_ok: function () {
    UI.restore_focus();
    $id('assistant_plotfunc').style.display = 'none';
    $id('assistant_plotfunc1var').style.display = 'none';
    var tmp = 'plotfunc(' + $id('plotfuncexpr').value + ',' + $id('plotfuncvarname').value + '===' + $id('plotfuncvarmin').value + '..' + $id('plotfuncvarmax').value;
    var tmp1 = $id('plotfuncvarstep').value;
    if (tmp1.length) tmp = tmp + ',xstep===' + tmp1;
    tmp += ')';
    UI.insert(UI.focused, tmp);
  },
  assistant_plotfunc2var_ok: function () {
    UI.restore_focus();
    $id('assistant_plotfunc').style.display = 'none';
    $id('assistant_plotfunc2var').style.display = 'none';
    var tmp = 'plotfunc(' + $id('plotfunc2expr').value + ',[' + $id('plotfunc2varx').value + ',' + $id('plotfunc2vary').value + ']';
    var tmp1 = $id('plotfunc2varxstep').value;
    if (tmp1.length) tmp = tmp + ',xstep===' + tmp1;
    tmp1 = $id('plotfunc2varystep').value;
    if (tmp1.length) tmp = tmp + ',ystep===' + tmp1;
    tmp += ')';
    UI.insert(UI.focused, tmp);
  },
  assistant_plotparam_ok: function () {
    UI.restore_focus();
    $id('assistant_plotparam').style.display = 'none';
    $id('assistant_plotparam1var').style.display = 'none';
    var tmp = 'plotparam([' + $id('plotparamexprx').value + ',' + $id('plotparamexpry').value + '],' + $id('plotparamvarname').value + '===' + $id('plotparamvarmin').value + '..' + $id('plotparamvarmax').value;
    var tmp1 = $id('plotparamvarstep').value;
    if (tmp1.length) tmp = tmp + ',tstep===' + tmp1;
    tmp += ',display===cap_flat_line)';
    UI.insert(UI.focused, tmp);
  },
  assistant_plotparam2var_ok: function () {
    UI.restore_focus();
    $id('assistant_plotparam').style.display = 'none';
    $id('assistant_plotparam2var').style.display = 'none';
    var tmp = 'plotparam([' + $id('plotparam2exprx').value + ',' + $id('plotparam2expry').value + ',' + $id('plotparam2exprz').value + '],[' + $id('plotparam2varx').value + ',' + $id('plotparam2vary').value + ']';
    var tmp1 = $id('plotparam2varxstep').value;
    if (tmp1.length) tmp = tmp + ',ustep===' + tmp1;
    tmp1 = $id('plotparam2varystep').value;
    if (tmp1.length) tmp = tmp + ',vstep===' + tmp1;
    tmp += ')';
    UI.insert(UI.focused, tmp);
  },
  assistant_plotimplicit_ok: function () {
    UI.restore_focus();
    $id('assistant_plotimplicit').style.display = 'none';
    var ctr = $id('plotimplicitlevel').value;
    var tmp = $id('plotimplicitexprf').value + ',[' + $id('plotimplicitvarx').value + ',' + $id('plotimplicitvary').value + ']';
    if (ctr.length) tmp = 'plotcontour(' + tmp + ',' + ctr; else tmp = 'plotimplicit(' + tmp;
    var tmp1 = $id('plotimplicitvarxstep').value;
    if (tmp1.length) tmp = tmp + ',xstep===' + tmp1;
    tmp1 = $id('plotimplicitvarystep').value;
    if (tmp1.length) tmp = tmp + ',ystep===' + tmp1;
    tmp += ')';
    UI.insert(UI.focused, tmp);
  },
  split_right:function(s,c='='){
    let p=s.search('=');
    while (p>=0 && p<s.length && s[p]=='=') ++p;
    if (p>=0 && p<=s.length)
      return s.substr(p,s.length-p);
    return s;
  },
  assistant_plotfield_ok: function () {
    UI.restore_focus();
    $id('assistant_plotfield').style.display = 'none';
    var ctr = $id('plotfieldinit').value;
    var tmp0 = 'gl_x='+UI.split_right($id('plotfieldvarx').value)+';gl_y='+UI.split_right($id('plotfieldvary').value)+';';
    var tmp = $id('plotfieldexprf').value + ',[' + $id('plotfieldvarx').value + ',' + $id('plotfieldvary').value + ']';
    if (ctr.length) tmp = tmp0+'plotfield(' + tmp + ',plotode===' + ctr; else tmp = tmp0+'plotfield(' + tmp;
    var tmp1 = $id('plotfieldvarxstep').value;
    if (tmp1.length) tmp = tmp + ',xstep===' + tmp1;
    tmp1 = $id('plotfieldvarystep').value;
    if (tmp1.length) tmp = tmp + ',ystep===' + tmp1;
    tmp += ')';
    UI.insert(UI.focused, tmp);
  },
  assistant_plotpolar_ok: function () {
    UI.restore_focus();
    var tmp = 'plotpolar(' + $id('plotpolarexpr').value + ',' + $id('plotpolarvarname').value + ',' + $id('plotpolarvarmin').value + ',' + $id('plotpolarvarmax').value;
    UI.insert(UI.focused, tmp);
    tmp = $id('plotpolarvarstep').value;
    if (tmp.length) UI.insert(UI.focused, ',tstep===' + tmp);
    UI.insert(UI.focused, ')');
    $id('assistant_plotpolar').style.display = 'none';
  },
  usecm: true,
  fixeddel: false,
  kbdshift: false,
  usemathjax: false,
  prettyprint: true,
  qa: false,
  histcount: 0,
  selection: '',
  langue: -1,
  calc: 1, // 1 KhiCAS, 2 Numworks, 3 TI Nspire CX
  xwaspy_shift: 33, // must be >32 for space encoding, and <=35 for a..z encoding
  //locked:false,
  canvas_w: 350,
  canvas_h: 200,
  canvas_lastx: 0,
  canvas_lasty: 0,
  canvas_pushed: false,
  touch_handler:function(event){
    var touches = event.changedTouches,
        first = touches[0];
    var s2=first.target.id;
    var is_3d= s2.length>5 && s2.substr(0,4)=='gl3d';
    var n3d='';
    if (is_3d)
      n3d = s2.substr(5, s2.length - 5);
    event.preventDefault();
    if (event.type=="touchstart"){
      UI.canvas_pushed=true;
      UI.canvas_lastx=first.clientX; UI.canvas_lasty=first.clientY;
    }
    if (event.type=="touchmove"){
      UI.canvas_mousemove(first, n3d);
    }
    if (event.type=="touchend"){
      UI.canvas_pushed=false;
    }
  },
  rendercomment: function (t) { // replace \n by <br>, for $$ call caseval on mathml or latex(quote())
    var res = '', i = 0, l = t.length, cas, indollar = false;
    // console.log(t,l);
    for (; i < l; ++i) {
      var ch = t.charAt(i), nxt = 0;
      if (i < l - 1) nxt = t.charAt(i + 1);
      if (ch == '$' && nxt == '$') {
        if (!UI.usemathjax) {
          var script = document.createElement("script");
          script.type = "text/javascript";
          script.src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML";
          document.getElementsByTagName("head")[0].appendChild(script);
          window.setTimeout(UI.setmathjax, 400);
        }
        return t;
      }
      // console.log(i,ch,nxt);
      if (ch == '\n') {
        if (!indollar) res += '<br>';
        continue;
      }
      if (ch == '\\' && nxt == '$') {
        if (indollar) cas += nxt; else res += nxt;
        continue;
      }
      if (ch != '$') {
        if (indollar) cas += ch; else res += ch;
        continue;
      }
      if (!indollar) {
        indollar = true;
        cas = '';
        continue;
      }
      indollar = false;
      if (UI.usemathjax)
        cas = 'latex(quote(' + cas + '))';
      else
        cas = 'mathml(quote(' + cas + '),1)';
      //console.log(cas);
      cas = UI.caseval_noautosimp(cas);
      //console.log(cas);
      if (cas.charAt(0) == '"' && cas.length > 2)
        cas = cas.substr(1, cas.length - 2);
      if (UI.usemathjax)
	cas='$$'+cas+'$$';
      res += cas;
    }
    return res;
  },
  add_math: function (field) {
    UI.insert(field, "\$ \$");
    UI.moveCaret(field, -2);
  },
  add_strong: function (field) {
    UI.insert(field, "<strong></strong>");
    UI.moveCaret(field, -9);
  },
  add_em: function (field) {
    UI.insert(field, "<em></em>");
    UI.moveCaret(field, -5);
  },
  add_red: function (field) {
    UI.insert(field, "<font color='red'></font>");
    UI.moveCaret(field, -7);
  },
  add_green: function (field) {
    UI.insert(field, "<font color='green'></font>");
    UI.moveCaret(field, -7);
  },
  add_tt: function (field) {
    UI.insert(field, "<tt></tt>");
    UI.moveCaret(field, -5);
  },
  add_href: function (field) {
    UI.insert(field, "<a href=\"http://\" target=\"_blank\"></a>");
    UI.moveCaret(field, -22);
  },
  add_img: function (field) {
    UI.insert(field, "<img width=\"32\" height=\"32\" src=\"\">");
    UI.moveCaret(field, -2);
  },
  add_list: function (field) {
    UI.insert(field, "<ul>\n<li>\n<li>\n<li>\n</ul>");
    UI.moveCaretUpDown(field, -3);
  },
  add_enum: function (field) {
    UI.insert(field, "<ol>\n<li>\n<li>\n<li>\n</ol>");
    UI.moveCaretUpDown(field, -3);
  },
  add_h1: function (field) {
    UI.insert(field, "<h1></h1>");
    UI.moveCaret(field, -5);
  },
  add_h2: function (field) {
    UI.insert(field, "<h2></h2>");
    UI.moveCaret(field, -5);
  },
  ckenter_comment: function (event, field) {
    var key = event.keyCode;
    //console.log('comment',event.keyCode,event.shiftKey,key);
    if (key==27){ event.preventDefault(); field.blur(); return;}
    if (key == 13) UI.resizetextarea(field);
    var skipline = key != 13 || !event.ctrlKey;
    if (skipline && key == 13 && !event.shiftKey) {
      // check enter pressed at beginning of field
      var pos = field.selectionStart;
      var end = field.selectionEnd;
      if (pos == 0 && end == 0)
        skipline = false;
      // check enter pressed at end after a blank line
      var text = field.value;
      if (pos == text.length && end == pos) {
        if (pos >= 1 && text[pos - 1] == '\n')
          skipline = false;
      }
    }
    if (skipline) return true;
    UI.editcomment_end(field.nextSibling, true);
    return false;
  },
  addcomment: function (text) {
    var s = '';
    var t = text; // text.length >= 3 ? text.substr(3, text.length - 3) : text;
    s += '<td colspan="2"><span style="display:none"><textarea title="';
    s += UI.langue == -1 ? 'Ctrl-Enter ou Ok pour valider ce commentaire' : 'Ctrl-Enter or Ok will update this comment';
    s += '" onfocus="UI.focused=this;UI.ckswitch(\'comment\',0)" onblur="UI.editcomment_end(this.nextSibling,false)" onkeydown="UI.ckenter_comment(event,this)" row="5" cols="60">' + t + '</textarea>';
    s += '<button class="bouton" onmousedown="event.preventDefault()" title="';
    s += UI.langue == -1 ? 'Valide ce commentaire' : 'Update comment';
    s += '" onclick="UI.editcomment_end(this,true)">Ok</button>';
    s += '<br><button onmousedown="event.preventDefault()" title="';
    s += UI.langue == -1 ? 'Ajoute $ $ pour ins&eacute;rer des maths' : 'Add $ $ to insert maths';
    s += '" class="bouton" onclick="UI.add_math(UI.focused)">math</button>';
    s += '<button onmousedown="event.preventDefault()" title="'
    s += UI.langue == -1 ? 'Passe en police fixe' : 'Fixed size police';
    s += '" class="bouton" onclick="UI.add_tt(UI.focused)"><tt>Abc</tt></button>';
    s += '<button title="';
    s += UI.langue == -1 ? 'Passe en gras' : 'Boldface';
    s += '" onmousedown="event.preventDefault()"  class="bouton" onclick="UI.add_strong(UI.focused)"><strong>Abc</strong></button>';
    s += '<button title="';
    s += UI.langue == -1 ? 'Passe en italique' : 'italics';
    s += '" onmousedown="event.preventDefault()" class="bouton" onclick="UI.add_em(UI.focused)"><em>Abc</em></button>';
    s += '<button title="Red/Rouge" onmousedown="event.preventDefault()" class="bouton" onclick="UI.add_red(UI.focused)"><font color="red">Abc</font></button>';
    s += '<button title="Green/Vert" onmousedown="event.preventDefault()" class="bouton" onclick="UI.add_green(UI.focused)"><font color="green">Abc</font></button>';
    s += '<button onmousedown="event.preventDefault()" title="';
    s += UI.langue == -1 ? 'Ins&egrave;re une liste' : 'List insertion';
    s += '" class="bouton" onclick="UI.add_list(UI.focused)">list</button>';
    s += '<button onmousedown="event.preventDefault()" title="';
    s += UI.langue == -1 ? 'Ins&egrave;re une listenum&eacute;rot&eacute;e' : 'Numbered list insertion';
    s += '" class="bouton" onclick="UI.add_enum(UI.focused)">num</button>';
    s += '<button onmousedown="event.preventDefault()" title="';
    s += UI.langue == -1 ? 'Ins&egrave;re un lien' : 'Link insertion';
    s += '" class="bouton" onclick="UI.add_href(UI.focused)">link</button>';
    s += '<button onmousedown="event.preventDefault()" title="';
    s += UI.langue == -1 ? 'Ins&egrave;re une image' : 'Image insertion';
    s += '" class="bouton" onclick="UI.add_img(UI.focused)">img</button>';
    s += '<button style="display:none" onmousedown="event.preventDefault()" class="bouton" onclick="UI.add_h2(UI.focused)">sous</button>';
    s += UI.langue == -1 ? '<button onmousedown="event.preventDefault()" title="Ins&egrave;re un titre" class="bouton" onclick="UI.add_h1(UI.focused)">titre</button>' : '<button onmousedown="event.preventDefault()" title="Insert title" class="bouton" onclick="UI.add_h1(UI.focused)">title</button>';
    s += '</span><span onclick="UI.editcomment1(this);">' + UI.rendercomment(t) + '</span></td>';
    s += '<td>';
    return s;
  },
  editcomment1: function (field) {
    //console.log('editcomment1');
    if (UI.readonly) return; 
    UI.tableaucm=1;
    var prev = field.parentNode, prevff = prev.firstChild.firstChild;
    //console.log(prevff);
    var pos = prevff.value.search('<a href="');
    if (pos >= 0 && pos < prevff.value.length)
      return;
    field.style.display = 'none';
    prev.firstChild.style.display = 'inline';
    // console.log(prev.firstChild.firstChild.innerHTML);
    UI.focused = prevff;
    UI.resizetextarea(UI.focused);
    UI.focused.focus();
    //console.log(prev);
    //prev.nextSibling.style.display = 'none';
  },
  editcomment_end: function (field, comment) {
    UI.tableaucm=0;
    var prev = field.previousSibling;
    var s = prev.value;
    if (s=='') s='#';
    var par = field.parentNode;
    //console.log(prev);
    par.style.display = 'none';
    par.nextSibling.style.display = 'inline';
    if (comment || par.nextSibling.innerHTML=='') {
      par.nextSibling.innerHTML = UI.rendercomment(s);
      if (UI.usemathjax)
        MathJax.Hub.Queue(["Typeset", MathJax.Hub, par.nextSibling]);
    }
    UI.link(0);
    UI.focused = 0;
  },
  giac_renderer: function (text) {
    //console.log('giac_renderer_1',text);
    var gr = Module.cwrap('_ZN4giac13giac_rendererEPKc', 'number', ['string']);
    gr(text);
    //console.log('giac_renderer_2',text);
    var keyboardListeningElement = Module['keyboardListeningElement'] || document;
    keyboardListeningElement.removeEventListener("keydown", SDL.receiveEvent);
    keyboardListeningElement.removeEventListener("keyup", SDL.receiveEvent);
    keyboardListeningElement.removeEventListener("keypress", SDL.receiveEvent);
  },
  canvas_mousemove: function (event, no) {
    if (UI.canvas_pushed) {
      // Module.print(event.clientX);
      if (UI.canvas_lastx != event.clientX) {
        if (event.clientX > UI.canvas_lastx)
          UI.giac_renderer('r' + no);
        else
          UI.giac_renderer('l' + no);
        UI.canvas_lastx = event.clientX;
      }
      if (UI.canvas_lasty != event.clientY) {
        if (event.clientY > UI.canvas_lasty)
          UI.giac_renderer('d' + no);
        else
          UI.giac_renderer('u' + no);
        UI.canvas_lasty = event.clientY;
      }
    }
  },
  initconfigstring: '',
  python_mode: 0,
  python_indent: 4,
  warnpy: false, // set to false if you do not want Python compat warning
  xtn: 'x', // var name, depends on last app
  listCookies: function () { // list cookies with name begin == ' xcas__'
    var aString = '';
    if (window.localStorage) {
      for (var i = 0, len = localStorage.length; i < len; i++) {
        var tmp = localStorage.key(i);
        //console.log(tmp);
        if (tmp.substr(0, 6) == 'xcas__') {
          var tmpname = tmp.substr(6, tmp.length - 6);
          aString += "<button onclick=\"UI.restorefrom('" + tmp.substr(0, pos) + "')\">" + tmpname + "</button>\n";
        }
      }
    }
    var theCookies = document.cookie.split(';');
    for (var i = 0; i < theCookies.length; i++) {
      // console.log(i,theCookies[i].substr(0,7));
      var tmp = theCookies[i];
      var pos = tmp.search('=');
      if (pos > 7 && tmp.substr(0, 7) == ' xcas__') {
        var tmpname = tmp.substr(7, pos - 7);
        aString += "<button onclick=\"UI.restorefrom('" + tmp.substr(1, pos - 1) + "')\">" + tmpname + "</button>\n";
      }
    }
    aString += "<button onclick=$id('loadfile_cookie').innerHTML=''>X</button>\n"
    //console.log(aString);
    return aString;
  },
  detectmob: function () {
    if (navigator.userAgent.match(/Android/i)
        || navigator.userAgent.match(/webOS/i)
        || navigator.userAgent.match(/iPhone/i)
        || navigator.userAgent.match(/iPad/i)
        || navigator.userAgent.match(/iPod/i)
        || navigator.userAgent.match(/BlackBerry/i)
        || navigator.userAgent.match(/Windows Phone/i)
    ) return true;
    else
      return false;
  },
  browser_type: function () {
    var isOpera = !!window.opera || navigator.userAgent.indexOf(' OPR/') >= 0;
    var isFirefox = typeof InstallTrigger !== 'undefined';   // Firefox 1.0+
    var isSafari = Object.prototype.toString.call(window.HTMLElement).indexOf('Constructor') > 0;
    var isChrome = !!window.chrome && !isOpera;              // Chrome 1+
    var isIE = /*@cc_on!@*/false || !!document.documentMode; // At least IE6
    if (isFirefox) return 1;
    if (isSafari) return 2;
    if (isChrome) return 3;
    if (isIE) return 4;
    if (isOpera) return 5;
    return 0;
  },
  is_touch_device: function () {
    return (('ontouchstart' in window)
        || (navigator.MaxTouchPoints > 0)
        || (navigator.msMaxTouchPoints > 0));
  },
  cmtotextarea:function(cm){
    var ta=cm.getTextArea();
    cm.toTextArea();
    UI.resizetextarea(ta);
    //console.log('cmtotextarea',ta.parentNode);
    ta.parentNode.style.top=UI.mobile_savetop;
    ta.parentNode.style.left=UI.mobile_saveleft;
    ta.parentNode.firstChild.style.display='inline';
    for (let nxt=ta.nextSibling;nxt;nxt=nxt.nextSibling){
      nxt.style.display='inherit';
    }
  },
  textareatocm:function(divText){
    UI.focused=divText;
    let prev=divText.previousSibling;
    if (prev && prev.tagName=='SPAN'){
      //console.log(prev.innerText);
      UI.ckswitch(prev.innerText,0);
    }
    for (let nxt=divText.nextSibling;nxt;nxt=nxt.nextSibling){
      nxt.style.display='none';
    }
    $id('assistant_keys').style.display='block';
    if (!UI.detectmob()){
      while (parseInt(divText.parentNode.style.left)>=window.innerWidth-BoardManager.scrollQuantity())
	BoardManager.right();
      UI.open_menu(false); // $id('assistant_keys').style.display='block';
      //console.log('textareatocm');console.trace();
    }
    UI.mobile_savetop=divText.parentNode.style.top;
    UI.mobile_saveleft=divText.parentNode.style.left;
    divText.parentNode.style.top='0px';
    if (UI.detectmob()){
      divText.parentNode.style.left='0px';
      // hide >
      divText.parentNode.firstChild.style.display='none';
    }
    if (!UI.usecm)
      return;
    var cmentree = CodeMirror.fromTextArea(divText, {
      matchBrackets: true,
      lineNumbers: true,
      lineWrapping: true,
      viewportMargin: Infinity,
    });
    UI.focused=cmentree;
    // cmentree.style.color = "rgba(0.2, 0.2, 0.2, 0.9)";
    // cmentree.style.backgroundColor = "white";
    UI.setoption_mode(cmentree);
    //console.log(entree.type);
    cmentree.options.indentUnit = UI.python_mode ? UI.python_indent : 2;
    UI.codemirror_setoptions(cmentree);
    cmentree.setValue(divText.value);
    if (UI.detectmob()){
      cmentree.setSize(0.9*window.innerWidth,0.23*window.innerHeight);
      let n=UI.count_newline(divText.value),f=18;
      if (n>=4 || n==5) f=17; 
      if (n>=6 && n<9) f=16;
      if (n>=9 && n<12) f=15;
      if (n>=12) f=14;
      UI.changefontsize(cmentree, f);
    }
    else {
      UI.changefontsize(cmentree, 20);
    }
    cmentree.focus(); UI.scroll_controls();
    return cmentree;
  },
  moveCaret: function (field, charCount) {
    if (field.type != "textarea" && field.type != "text") {
      var pos = field.getCursor();
      pos.ch = pos.ch + charCount;
      field.setCursor(pos);
      field.refresh();
      // UI.show_curseur();
      return;
    }
    var pos = field.selectionStart;
    pos = pos + charCount;
    if (pos < 0) pos = 0;
    if (pos > field.value.length) pos = field.value.length;
    field.setSelectionRange(pos, pos);
  },
  backspace: function (field) {
    //if (UI.focusaftereval) field.focus();
    if (field.type != "textarea" && field.type != "text") {
      var start = field.getCursor('from');
      var end = field.getCursor('to');
      if (end.line != start.line || end.ch != start.ch)
        field.replaceSelection('');
      else {
        var c = start.ch;
        var l = start.line;
        if (start.ch == 0 && start.line == 0){
	  return;
	}
        if (c > 0) {
          var s = field.getRange({line: l, ch: 0}, end), i;
          for (i = 0; i < s.length; i++) {
            if (s.charAt(i) != ' ') break;
          }
          //console.log(i,s.length,c);
          if (i == s.length && c >= 2) {
            var l1 = l - 1;
            for (; l1 >= 0; --l1) {
              s = field.getLine(l1);
              for (i = 0; i < s.length && i < c; i++) {
                if (s[i] != ' ') break;
              }
              if (i != s.length && i < c) break;
            }
            if (l1 >= 0) c = i; else c -= 2;
          }
          else c--;
          field.replaceRange('', {line: l, ch: c}, end);
        }
        else {
          l--;
          var s = field.getRange({line: l, ch: c}, end);
          field.replaceRange('', {line: l, ch: s.length - 1}, end);
        }
      }
      var t = field.getTextArea();
      t.value = field.getValue();
    } else {
      var pos = field.selectionStart;
      var pos2 = field.selectionEnd;
      var s = field.value;
      if (pos < pos2) {
        field.value = s.substring(0, pos) + s.substring(pos2, s.length);
        if (pos < 0) pos = 0;
        if (pos > field.value.length) pos = field.value.length;
        field.setSelectionRange(pos, pos);
        UI.resizetextarea(field);
        return;
      }
      if (pos > 0) {
        field.value = s.substring(0, pos - 1) + s.substring(pos, s.length);
        pos--;
        if (pos < 0) pos = 0;
        if (pos > field.value.length) pos = field.value.length;
        field.setSelectionRange(pos, pos);
        UI.resizetextarea(field);
      }
    }
  },
  insert_focused:function(value,close=2){
    //console.log(value);
    if (close==2 ||
	(close==1 && UI.detectmob())
       ){
      $id('assistant_cas').style.display='none';
      UI.close_menus();
    }
    UI.insert(UI.focused,value);
    UI.scroll_controls();
    if (UI.focused!=0) UI.focused.focus();
  },
  insert: function (field, value) {
    var myValue = value.replace(/&quote;/g, '\"');
    //console.log('insert',field);
    if (UI.focusaftereval){ field.focus(); UI.scroll_controls(); }
    if (field.type != "textarea" && field.type != "text") {
      if (field.type == undefined && field.firstChild != undefined) {
        // console.log(field.innerHTML);
        UI.insert(field.firstChild, value);
        return;
      }
      var start = field.getCursor('from');
      var end = field.getCursor('to');
      if (end.line != start.line || end.ch != start.ch || myValue.length < 3)
        field.replaceSelection(myValue);
      else {
        // detect the same command not selected just before
        var parpos = 0;
        for (; parpos < myValue.length; ++parpos) {
          if (myValue[parpos] == '(') break;
        }
        if (parpos < myValue.length) {
          var S = {line: start.line, ch: start.ch};
          S.ch -= (parpos + 1);
          var deb = S.ch == -1;
          if (deb) S.ch = 0;
          if (S.ch >= 0) {
            var avant = field.getRange(S, end);
            //console.log('avant avant',avant,S.ch,deb);
            var tst = avant.length - 1;
            if (!deb && avant[tst] != '(') {
              avant = avant.substr(1, avant.length - 1);
              //console.log('avant avant',avant);
              --tst;
            }
            else {
              if (avant[tst] == '(')
                ++parpos;
            }
            //console.log('avant apres',avant);
            for (; tst >= 0; --tst) {
              if (myValue[tst] != avant[tst]) break;
            }
            //console.log('avant',myValue);
            if (tst < 0)
              myValue = myValue.substr(parpos, myValue.length - parpos);
            //console.log('apres',myValue);
          }
        }
        field.replaceSelection(myValue);
      }
      field.execCommand("indentAuto");
      var t = field.getTextArea();
      t.value = field.getValue();
    }
    else {
      var pos = field.selectionStart;
      pos = pos + myValue.length;
      //IE support
      if (document.selection) {
        if (UI.focusaftereval){ field.focus(); UI.scroll_controls();}
        var sel = document.selection.createRange();
        sel.text = myValue;
      }
      //MOZILLA and others
      else {
        var startPos = field.selectionStart;
        var endPos = field.selectionEnd;
        if (field.selectionStart || field.selectionStart == '0') {
          field.value = field.value.substring(0, startPos)
              + myValue
              + field.value.substring(endPos, field.value.length);
        } else {
          field.value += myValue;
        }
      }
      field.setSelectionRange(pos, pos);
      UI.resizetextarea(field);
    }
    // UI.show_curseur();
  },  
  codemirror_setoptions: function (cmentree) {
    UI.setoption_mode(cmentree);
    cmentree.on("blur", function (cm) {
      if (UI.savefocused || UI.hints || $id('online_help').style.display!='none') return;
      // console.log('blur',cm);
      UI.tableaucm=0; UI.clear_focus(); 
      if (cm.getSelection().length > 0) {
        UI.selection = cm.getSelection();
      }
      var ta=cm.getTextArea();
      var v=ta.value;
      UI.cmtotextarea(cm);
      ta.innerHTML=ta.value=v;
    });
    cmentree.on("focus", function (cm) {
      // console.log('focus',cm);
      UI.tableaucm=1;
      //UI.set_config_width();
    });
    cmentree.setOption("extraKeys", {
      Enter: function (cm) { // guess if newline evaluates or adds a newline
        var start = cm.getCursor('from');
        var end = cm.getCursor('to');
        var tst = cm.lineCount() > 1;
        if (!tst) { // if the line begins by function/fonction or def/for/if/while/si/tantque/pour
          var txt = cm.getLine(end.line);
          //console.log(txt);
          var l = txt.length, i, j;
          for (i = 0; i < l; i++) {
            if (txt[i] != ' ') break;
          }
          for (j = i; j < l; j++) {
            if (txt[j] == ' ') break;
          }
          txt = txt.substr(i, j - i);
          //console.log(txt);
          if (txt == "for" || txt == "while" || txt == "if" || txt == "pour" || txt == "tantque" || txt == "si" || txt == "def" || txt == "fonction" || txt == "function")
            tst = true;
        }
        if (tst && (end.line != start.line || end.ch != start.ch || ((start.line > 0 || start.ch > 0) && UI.not_empty(cm.getLine(end.line)))))
          UI.insert(cm, '\n');
        else {
	  UI.tableaucm=-1;
	  var ta=cm.getTextArea();
	  ta.innerHTML=ta.value=cm.getValue();
	  var taille=(parseInt(ta.style.fontSize)*(UI.count_newline(ta.value)+1)+6);
	  // console.log('taille',taille);
	  ta.style.height=taille+'px';
	  UI.cmtotextarea(cm);
          //UI.set_editline(cmentree,false);
          UI.tableaueval(ta);
	}
      },
      "Ctrl-Enter": function (cm) {
        //UI.set_editline(cmentree,false);
	UI.tableaucm=-1;
	var ta=cm.getTextArea();
	var v=cm.getValue();
	var taille=(parseInt(ta.style.fontSize)*(UI.count_newline(ta.value)+1)+6);
	//console.log('taille',taille);
	ta.style.height=taille+'px';
	UI.cmtotextarea(cm);
        ta.innerHTML=ta.value=v;
	UI.tableaueval(ta);
      },
      Backspace: function (cm) {
        UI.backspace(cm);
      },
      F1: function (cm) {
        UI.completion(cm);
      },
      "Ctrl-Space": function (cm) {
        UI.completion(cm);
      },
      Tab: function (cm) {
        UI.indent_or_complete(cm);
      },
      Esc: function (cm){
	UI.tableaucm=-1;
	var ta=cm.getTextArea();
	cm.setValue(ta.value);
	UI.cmtotextarea(cm);
      },
    });
  },
  setoption_mode: function (cmentree) {
    if (!UI.usecm) return;
    if (UI.python_mode) {
      //console.log('Python mode');
      if (UI.micropy)
	cmentree.setOption("mode", "micropy");
      else
	cmentree.setOption("mode", "python");
    }
    else {
      //console.log('Xcas mode');
      cmentree.setOption("mode", "simplemode");
    }
  },
  add_autosimplify: function (value) {
    var n = value.search(';');
    if (value.length == 0 || value[0] == '@')
      return value;
    if (n < 0 || n >= value.length) {
      var n = value.search('//');
      if (n < 0 || n >= value.length) {
        if (UI.python_mode)
          return 'add_autosimplify(@@' + value + ')';
        return 'add_autosimplify(' + value + ')';
      }
    }
    return value;
  },
  // emscripten micropython interface
  mp_init:function(taille){
    var init = Module.cwrap('mp_js_init', 'null', ['number']);
    UI.micropy_initialized=1;
    return init(taille);
  },
  mp_str:function(s){
    var ev = Module.cwrap('mp_js_do_str', 'number', ['string']);
    return ev(s);
  },
  micropy:0,
  micropy_initialized:0,
  micropy_heap:4194304,
  python_output:"",
  clean: function (text, quote) {
    var cmd = text;
    if (quote) cmd = cmd.replace(/\'/g, '\\\'');
    cmd = cmd.replace(/>/g, '&gt;');
    cmd = cmd.replace(/</g, '&lt;');
    if (quote) cmd = cmd.replace(/\"\"/g, '&quote;');
    return cmd;
  },
  clean_for_html: function(text){
    text = text.replace(/&/g, "&amp;");
    text = text.replace(/</g, "&lt;");
    text = text.replace(/>/g, "&gt;");
    text = text.replace(/\n/g, '<br>');
    return text;
  },    
  add_python_output:function(s){
    UI.python_output += s;
    //console.log(s);//console.log(UI.python_output);
  },
  set_config_width: function(){
  },
  sleep: function (miliseconds) {
    var currentTime = new Date().getTime();
    while (currentTime + miliseconds >= new Date().getTime()) {
    }
  },
  quickjs:function(text){
    while (text.length>0){
      var ch=text.substr(text.length-1,1);
      if (ch!=' ')
	break;
      text=text.substr(0,text.length-1);
    }
    if (text=='xcas' || text=='xcas '){
      UI.micropy=0; UI.python_mode=0; 
      var form = $id('config');
      form.python_xor.checked = false;
      form.python_mode.checked = false;
      form.js_mode.checked=false;
      UI.set_settings();
      return UI.caseval('python_compat(0)');
    }
    if (text=='.'){ // show turtle
      let s=UI.caseval('avance(0)');
      //console.log(s);
      return s;
    }
    if (text==','){ // show (matplotl)
      Module.print('>>> show()');
      let s=UI.caseval('show()');
      return s;
    }
    if (text==';'){
      let s=UI.caseval('show_pixels()');
      return s;
    }
    if (text.length>=2 && text[0]=='@'){
      if (text[1]=='@')
	return eval(text.substr(2,text.length-2));
      text=text.substr(1,text.length-1);
    }
    else text='"use math";'+text;
    let ev=Module.cwrap('quickjs_ck_eval', 'string', ['string']);
    return ev(text);
  },  
  mp_eval:function(text){
    while (text.length>0){
      var ch=text.substr(text.length-1,1);
      if (ch!=' ')
	break;
      text=text.substr(0,text.length-1);
    }
    if (text=='xcas' || text=='xcas '){
      UI.micropy=0; UI.python_mode=0; 
      //var form = $id('config');
      //form.python_xor.checked = false;
      //form.python_mode.checked = true;
      //UI.set_settings();
      return UI.caseval('python_compat(1)');
    }
    if (text=='.'){ // show turtle
      var s=UI.caseval('avance(0)');
      //console.log(s);
      return s;
    }
    if (text==','){ // show (matplotl)
      Module.print('>>> show()');
      var s=UI.caseval('show()');
      return s;
    }
    if (text==';'){
      var s=UI.caseval('show_pixels()');
      return s;
    }
    if (!UI.micropy_initialized){
      UI.mp_init(UI.micropy_heap);
      console.log('mp_init done');
      try {
	UI.mp_eval("1");
      }
      catch(err) { console.log('mp_init eval(1)',err);}
    }
    UI.python_output='';
    /*
    var pos=text.search('=');
    if (pos<0){
      pos=text.search('print');
      if (pos<0)
	text='print('+text+')';
    }
    */
    console.log('mp_eval',text);
    //Module.print('>>> '+text);
    UI.mp_str(text);
    console.log('mp_evaled',UI.python_output);
    if (UI.python_output==''){
      return '"Done"';
    }
    if (UI.python_output.substr(UI.python_output.length-1,1)=='\n')
      UI.python_output=UI.python_output.substr(0,UI.python_output.length-1);
    if (UI.python_output.length>4 && UI.python_output.substr(0,5)=='"<svg')
      return UI.caseval('show()');
    return UI.python_output; // '"'+UI.python_output+'"';
  },
  handle_shortcuts:function(text){
    while (text.length>0){
      var ch=text.substr(text.length-1,1);
      if (ch!=' '){
	if (ch==':') return text.substr(0,text.length-1)+';show_pixels()';
	break;
      }
      text=text.substr(0,text.length-1);
    }
    if (text=='.') return 'avance(0)';
    if (text==',') return 'show()';
    if (text==';') return 'show_pixels()';
    if (text=='python' || text=='python '){
      UI.micropy=1; UI.python_mode = 4;
      //UI.set_settings();
      //var form = $id('config');
      //form.python_xor.checked = true;
      //form.python_mode.checked = true;
      console.log('now eval with micropython');
      return 'python_compat(4)';
    }
    if (text=='javascript' || text=='javascript '){
      UI.micropy=-1; UI.python_mode = 0;
      console.log('now eval with javascript');
      return 'python_compat(0)';
    }
    return text;
  },
  caseval: function (text) {
    if (!UI.ready) return ' Not ready ';
    var docaseval = Module.cwrap('caseval', 'string', ['string']);
    if (!UI.colorinit){
      console.log('default color',docaseval('color(white)'));
      UI.colorinit=true;
    }
    var value = UI.handle_shortcuts(text);
    value = value.replace(/%22/g, '\"');
    value = UI.add_autosimplify(value);
    var s, err;
    // console.log('caseval',value);
    try {
      s = docaseval(value);
    } catch (err) {
      console.log('caseval error ',err);
    }
    // Module.print(text+ ' '+s);
    return s;
  },
  caseval_noautosimp: function (text) {
    if (!UI.ready) return ' Clic_on_Exec ';
    //console.log(text);
    var docaseval = Module.cwrap('caseval', 'string', ['string']);
    var value = UI.handle_shortcuts(text);
    value = value.replace(/%22/g, '\"');
    var s, err;
    try {
      s = docaseval(value);
    } catch (err) {
      console.log(err);
    }
    return s;
  },
  ckswitch:function(s1,al=1){
    //console.log('ckswitch',s1); console.trace();
    if (s1.length>1 && s1[s1.length-1]==' ')
      s1=s1.substr(0,s1.length-1);
    if (s1=='khicas' || s1=='cas' || s1=='>') s1='xcas';
    if (s1=='js' || s1=='>>') s1='javascript';
    if (s1=='micropy' || s1=='py' || s1=='>>>') s1='python';
    if (s1=='javascript'){
      UI.magnet_evaluator=s1;
      UI.createCookie('xcas_magnet_button',s1);
      if (al)
	alert('Javascript eval');
      $id('eval_js').style.fontWeight="900";
      $id('eval_xcas').style.fontWeight="normal";
      $id('eval_mp').style.fontWeight="normal";
      $id('eval_comment').style.fontWeight="normal";
      $id('addmagnetbutton').style.backgroundColor="orange";
      $id('addmagnetbutton').style.color="black";
      $id('addmagnetbutton').innerHTML="JS";
      $id('addmagnetbutton2').style.backgroundColor="orange";
      $id('addmagnetbutton2').style.color="black";
      $id('addmagnetbutton2').innerHTML="JavaScript";
      UI.micropy=-1;UI.python_mode = 0;
      return UI.caseval('python_compat(0)');
    }
    if (s1=='xcas'){
      UI.magnet_evaluator=s1;
      UI.createCookie('xcas_magnet_button',s1);
      if (al)
	alert('Xcas eval');      
      $id('addmagnetbutton').style.backgroundColor="cyan";
      $id('addmagnetbutton').style.color="black";
      $id('addmagnetbutton').innerHTML="Xc";
      $id('addmagnetbutton2').style.backgroundColor="cyan";
      $id('addmagnetbutton2').style.color="black";
      $id('addmagnetbutton2').innerHTML="Xcas";
      $id('eval_xcas').style.fontWeight="900";
      $id('eval_js').style.fontWeight="normal";
      $id('eval_mp').style.fontWeight="normal";
      $id('eval_comment').style.fontWeight="normal";
      UI.micropy=0;UI.python_mode = 0;
      return UI.caseval('python_compat(0)');
    }
    if (s1=='python'){
      UI.magnet_evaluator=s1;
      UI.createCookie('xcas_magnet_button',s1);
      if (al)
	alert('Python eval');
      $id('addmagnetbutton').style.backgroundColor="yellow";
      $id('addmagnetbutton').style.color="black";
      $id('addmagnetbutton').innerHTML="Py";
      $id('addmagnetbutton2').style.backgroundColor="yellow";
      $id('addmagnetbutton2').style.color="black";
      $id('addmagnetbutton2').innerHTML="Python";
      $id('eval_xcas').style.fontWeight="normal";
      $id('eval_js').style.fontWeight="normal";
      $id('eval_mp').style.fontWeight="900";
      $id('eval_comment').style.fontWeight="normal";
      UI.micropy=1;UI.python_mode = 4;
      return UI.caseval('python_compat(4)');
    }
    if (s1=='comment'){
      UI.magnet_evaluator=s1;
      UI.createCookie('xcas_magnet_button',s1);
      if (al)
	alert('Commentaires texte/$maths$');
      $id('addmagnetbutton').style.backgroundColor="rgb(108, 108, 255)";
      $id('addmagnetbutton').style.color="white";
      $id('addmagnetbutton').innerHTML="Txt";
      $id('addmagnetbutton2').style.color="white";
      $id('addmagnetbutton2').style.backgroundColor="rgb(108, 108, 255)";
      $id('addmagnetbutton2').innerHTML="Txt";
      $id('eval_xcas').style.fontWeight="normal";
      $id('eval_js').style.fontWeight="normal";
      $id('eval_mp').style.fontWeight="normal";
      $id('eval_comment').style.fontWeight="900";
      return -2;
    }
    return -2;
  },
  addmagnet:function(x,y,evaluator='',content='',evalmode=0){
    //console.log('addmagnet',x,y,evaluator,content,evalmode);
    // eval==0 will switch to codemirror and select the text, ==1 eval the magnet, ==-1 nothing
    //console.log('addmagnet',x,y,evaluator);
    if (evaluator=='micropy')
      evaluator='py';
    if (!UI.readonly)
      UI.scroll_controls();
    if (UI.usecm)
      MagnetManager.addMagnetText(x, y,'cm',evaluator,content,evalmode);
    else
      MagnetManager.addMagnetText(x, y,'textarea',evaluator,content,evalmode);
  },
  magnet_evaluator:'comment',
  textmagnet_onkeydown: function(e,divText)  {
    // UI.focused=divText;
    let setFontSize = (size) => {
      divText.style.fontSize = size + "px";
      for (let o of divText.children) {
	o.style.fontSize = size + "px";
      }
    }
    if (e.key == "Escape") {
      divText.blur();
      window.getSelection().removeAllRanges();
      /*if(divText.innerHTML == "")
	MagnetManager.remove(div);*/
    }
    if ( (e.ctrlKey || divText.tagName=='DIV') && e.key=="Enter"){
      e.stopPropagation();
      UI.tableaueval(divText);
      return;
    }
    if ((e.ctrlKey && e.key == "=") || (e.ctrlKey && e.key == "+")) { // Ctrl + +
      let size = parseInt(divText.style.fontSize);
      size++;
      setFontSize(size);
      e.preventDefault();
    }
    if (e.ctrlKey && e.key == "-") { // Ctrl + -
      let size = parseInt(divText.style.fontSize);
      if (size > 6) size--;
      setFontSize(size);
      e.preventDefault();
    }
    e.stopPropagation();
  },
  // load/save
  remove_extension: function(name){
    var s=name.length,i;
    for (i=s-1;i>=0;--i){
      if (name[i]=='.')
	break;
    }
    if (i>0)
      return name.substr(0,i);
    return name;
  },
  svg_counter: 0,
  savesvg: function (field) {
    var s = field.innerHTML;
    s = '<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n' + s;
    var blob = new Blob([s], {type: "text/plain;charset=utf-8"});
    filename = $id("name").value;
    ++UI.svg_counter;
    filename += UI.svg_counter + ".svg";
    saveAs(blob, filename);
  },
  tableaufixload:function(magnets){
    for (let i = 0; i < magnets.length; i++) {
      var cur=magnets[i];
      UI.remove_3d_canvas(cur);
      if (cur.tagName=='DIV' && cur.hasChildNodes()){
	if (cur.firstChild.type=='textarea'){
	  //console.log('fixload',cur.firstChild);
	  // make it convertible to codemirror
	  let divText=cur.firstChild;
	  if (UI.usecm) divText.addEventListener('focus',()=>{UI.textareatocm(divText);});
	  divText.addEventListener('blur',()=>{if (!UI.usecm){ divText.parentNode.style.top=UI.mobile_savetop;divText.parentNode.style.left=UI.mobile_saveleft;}});
	}
	if (cur.firstChild.nextSibling && cur.firstChild.nextSibling.type=='textarea'){
	  //console.log('fixload',cur.firstChild.nextSibling);
	  // make it convertible to codemirror
	  let divText=cur.firstChild.nextSibling;
	  if (UI.usecm) divText.addEventListener('focus',()=>{UI.textareatocm(divText);});
	  divText.addEventListener('blur',()=>{if (!UI.usecm){ divText.parentNode.style.top=UI.mobile_savetop;divText.parentNode.style.left=UI.mobile_saveleft;}});
	  divText.onpointerdown = (e) => {  e.stopPropagation(); }
	  divText.onpointermove = (e) => { e.stopPropagation(); }
	  divText.onpointerup = (e) => {  e.stopPropagation(); }
	  divText.onkeydown = (e) => { UI.textmagnet_onkeydown(e,divText); }
	}
      }
      UI.render_canvas(cur);
    }    
    var ua = window.navigator.userAgent;
    var old_ie = ua.indexOf('MSIE ');
    var new_ie = ua.indexOf('Trident/');
    if ((old_ie > -1) || (new_ie > -1) ) { UI.usecm=false;}
  },
  tableaucrashrestore:function(name){
    if (window.localStorage){
      container.scrollTo({ top: 0, left: 0 }); // behavior: 'smooth'
      let v=localStorage.getItem(name);
      LoadSave.loadJSON(JSON.parse(v));
      v=localStorage.getItem(name+'.stack');
      if (v){
	let hashParams=v.split('&');
	console.log(hashParams);
	UI.parse(hashParams,true);
	BoardManager.redo();
      }
    }
  },
  tableaulocalsave:function(name){
    UI.link(0);
    if (window.localStorage){
      // first save current tableau in localstorage, in case of a crash
      const x1 = container.scrollLeft;
      BoardManager.shift_magnets(x1); // move magnets to their position if container not scrolled
      let magnets = document.getElementById("magnets").innerHTML;
      BoardManager.shift_magnets(-x1); // restore position
      let canvasDataURL = document.getElementById("canvas").toDataURL();
      let obj = { magnets: magnets, canvasDataURL: canvasDataURL };
      localStorage.setItem(name, JSON.stringify(obj));
      localStorage.setItem(name+'.stack',UI.handwritings_string());
    }
  },
  decompact:function(ch){
    let c=ch.charCodeAt(0);
    if (ch>='0' && ch<='9')
      return c-48;
    if (ch>='A' && ch<='Z')
      return c-55;
    if (ch>='a' && ch<='z')
      return c-61;
  },
  _compact:function(n){
    if (n>=0 && n<10)
      return String.fromCharCode(n+48);
    if (n>=10 && n<36)
      return String.fromCharCode(n+55);
    if (n>=36 && n<62)
      return String.fromCharCode(n+61);
  },
  compact:function(n){
    // n must be in [0,62^3-1]
    let n1=n%62,n2=((n-n1)/62)%62,n3=(n-n1-n2*62)/(62*62);
    return UI._compact(n3)+UI._compact(n2)+UI._compact(n1);
  },
  tableauloadlink:function(chaine){
    var hashParams = chaine.split('&');
    // console.log(hashParams);
    if (hashParams.length==0) return;
    if (hashParams.length==1 && chaine[0]!='='){
      // old format handwritings_string + magnets will contain at least 2 fields exec mode and radian mode
      if (chaine.length<4 || chaine.substr(0,4)!='http'){
	var s=window.location.href;
	if (s.length>4 && s.substr(0,4)=='file'){
	  container.scrollTo({ top: 0, left: 0 });
	  //console.log(chaine);
	  var v=localStorage.getItem(chaine);
	  if (v!=null){
	    Menu.hide();
	    Welcome.hide();
	    LoadSave.loadJSON(JSON.parse(v));
	    return;
	  }
	  alert('Pour des raisons de securite, vous devez charger le fichier '+chaine+' depuis le menu Charger Parcourir');
	  Menu.show();
	  return;
	}
	s=s.substr(0,s.indexOf('#'));
	if (s.length>5 && s.substr(s.length-5,5)=='.html'){
	  for (var l=s.length-1;l>=0;--l){
	    if (s[l]=='/'){
	      s=s.substr(0,l+1);
	      break;
	    }
	  }
	}
	chaine=s+chaine;
      } // end if chaine is not http
      $id("welcome").style.display='hidden';
      welcome.hidden=true; Welcome.hide();
      console.log('parameter',chaine);
      UI.httpload(chaine);
      if (0 && UI.detectmob()){
	$id("mobile_leftright").style.display='inline';
	controls.hidden=true;
	//document.getElementById("viewport").setAttribute("content", "width=device-width, initial-scale=0.75, maximum-scale=1.0, minimum-scale=0.25, user-scalable=yes");
      }
      UI.set_readonly(true);
      //Menu.hide();
    }
    // new format
    container.scrollTo({ top: 0, left: 0 }); // behavior: 'smooth'
    UI.parse(hashParams,false);
    BoardManager.redo();
  },
  radian_mode:1,
  radian:function(r){
    UI.radian_mode=r;
    UI.createCookie('xcas_angle_radian',UI.radian_mode);
    $id('angle_radian_yes').style.fontWeight=r==1?"900":"normal";
    $id('angle_radian_no').style.fontWeight=r==1?"normal":"900";
    window.setTimeout(UI.caseval_noautosimp, 100, 'angle_radian:='+r);    
  },    
  allow_exec:true, // set to false to disable execution of magnets from link load
  allow:function(b){
    UI.allow_exec=b;
    UI.createCookie('xcas_magnet_allow_exec',UI.allow_exec);
    $id('allow_exec_yes').style.fontWeight=b?"900":"normal";
    $id('allow_exec_no').style.fontWeight=b?"normal":"900";
  },
  parse_util:function(s){
    if (s[0]>='0' && s[0]<='7'){
      // optimized encoding for free drawing: s[0] 0 to 7 is the chalk color
      // followed by 3 char and 3 char for x and y position
      let colors = UI.palette_colors;
      let color=colors[s.charCodeAt(0)-48];
      let x=(UI.decompact(s[1])*62+UI.decompact(s[2]))*62+UI.decompact(s[3]);
      let y=(UI.decompact(s[4])*62+UI.decompact(s[5]))*62+UI.decompact(s[6]);
      s='='+color+","+x+","+y+","+s.substr(7,s.length-7);
    }
    return s;
  },
  parse_int:function(s){
    let i=0,l=s.length,r=0;
    for (;i<l;++i){
      let c=s.charCodeAt(i);
      if (c<48 || c>57){
	console.log('invalid char in '+s+' position '+i);
	console.trace();
	return 0;
      }
      r=r*10+(c-48);
    }
    return r;
  },
  parse:function(hashParams,restorecrashed=false){
    let doexec=false;
    for (var i = 0; i < hashParams.length; i++) {
      let s=UI.parse_util(hashParams[i]);
      if (s == 'exec') {
        console.log(s);
        if (UI.allow_exec || !restorecrashed)
	  doexec = true;
        continue;
      }
      if (s.length>4 && s.substr(0,4)=='+///')
	s='comment='+s.substr(4,s.length-4);
      // s = s.replace(/___/g, '%');
      let p=s.split('=');
      // s = s.replace(/%3b/g, ';');
      if (p[0]=='width'){
	canvas.width=UI.parse_int(p[1]);
	continue;
      }
      if (p[0]=='height'){
	canvas.height=UI.parse_int(p[1]);
	continue;
      }
      //console.log(p);
      if (p[0]=='') {
	//UI.draw_encoded(p[1],document.getElementById("canvas").getContext("2d"));
	if (p[1]){
	  //console.log('restore in cancelstack',p[1]);
	  BoardManager.cancelStack.push(p[1]);
	  UI.set_readonly(true);
	}
	continue;
      }
      if (p[0]=='python'){
	let mode=UI.parse_int(p[1]);
	if (mode==-1){
	  UI.ckswitch('javascript',0);
	}
	else {
	  UI.python_mode=mode;
	  if (mode==4)
	    UI.ckswitch('python',0);
	}
      }
      if (p[0]=='blackboard'){
	if (palette.colors[0]!=p[1])
	  BlackVSWhiteBoard.switch();
      }
      if (p[0]=='textarea') p[0]='cas'; // old saves
      if (p[0]=='cas' || p[0]=='py' || p[0]=='micropy' || p[0]=='js' || p[0]=='comment' || p[0]=='handwriting' || p[0]=='svg' || p[0]=='img'){
	let ms = decodeURIComponent(p[1]); // console.log(p[1]);
	if (ms.length && ms[0]==',')
	  ms=ms.substr(1,ms.length-1);
	//console.log(ms);
	let pos=ms.search(',');
	let ms0=ms.substr(0,pos),msbgcol='',mscol='';
	// background color
	// rgb(,,) or rgba(,,,) color, skip to )
	if (pos>4 && ms0.substr(0,3)=='rgb'){
	  pos=ms.indexOf(')')+1; // end of rgba(
	}
	if (ms.length && ms[0]>'9'){
	  msbgcol=ms.substr(0,pos);
	  ms=ms.substr(pos+1,ms.length-pos-1);
	  pos=ms.search(',');
	  ms0=ms.substr(0,pos);
	}
	// color, rgb skip
	if (pos>4 && ms0.substr(0,3)=='rgb'){
	  pos=ms.indexOf(')')+1; // end of rgba(
	}
	if (ms.length && ms[0]>'9'){
	  mscol=ms.substr(0,pos);
	  ms=ms.substr(pos+1,ms.length-pos-1);
	  pos=ms.search(',');
	  ms0=ms.substr(0,pos);
	}
	let mx=UI.parse_int(ms0);
	ms=ms.substr(pos+1,ms.length-pos-1);
	pos=ms.search(',');
	let my=UI.parse_int(ms.substr(0,pos));
	ms=ms.substr(pos+1,ms.length-pos-1);
	//console.log('restore eval magnet',p[0],mx,my,ms);
	UI.addmagnet(mx,my,p[0],ms,doexec?1:-1);
	let curm=MagnetManager.currentMagnet;
	curm.style.left=mx+'px';
	if (mscol!='')
	  curm.style.color=mscol;
	if (msbgcol!=''){
	  curm.style.backgroundColor=msbgcol;
	  if (msbgcol=='white' && mscol=='')
	    curm.style.color='black';
	}
	continue;
      }
      if (p[0]=='exec' && p[1]=='1') {
        console.log(s);
        doexec = true;
        continue;
      }
      if (p[0]=='radian') {
        if (p[1]=='0') {
	  UI.radian(0);
          // form.angle_mode.checked = false;
          window.setTimeout(UI.caseval_noautosimp, 100, 'angle_radian:=0');
          // window.setTimeout(UI.set_settings, 300);
        }
        if (p[1]=='1') {
	  UI.radian(1);
          //form.angle_mode.checked = true;
        }
        continue;
      }
    } // end for
  },    
  /*
  do_popstate:1,
  popstate:function(event) {
    // The popstate event is fired each time when the current history entry changes.
    var r = UI.do_popstate==0?true:confirm("OK: Quit");
    if (r == true) {
      UI.do_popstate=0;
      // Call Back button programmatically as per user confirmation.
      history.back();
      // Uncomment below line to redirect to the previous page instead.
      // window.location = document.referrer // Note: IE11 is not supporting this.
    } else {
      // Stay on the current page.
      $id("controls").scrollIntoView();
      Menu.show();
    }
  }, 
  */
  httpload:function (theUrl){
    if (window.XMLHttpRequest)
    {// code for IE7+, Firefox, Chrome, Opera, Safari
      xmlhttp=new XMLHttpRequest();
    }
    else
    {// code for IE6, IE5
      xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
    }
    xmlhttp.onreadystatechange=function()
    {
      if (xmlhttp.readyState==4 && xmlhttp.status==200)
      {
        var s=xmlhttp.responseText;
	console.log(s.substr(0,100));
	LoadSave.loadJSON(JSON.parse(s));
      }
    }
    xmlhttp.open("GET", theUrl, false );
    xmlhttp.send();    
  },
  // eval
  mobile_savetop:'0px',
  mobile_saveleft:'0px',
  tableaueval:function(cur){
    let divText=cur;
    if (cur.tagName=='DIV' && cur.hasChildNodes()){
      if (cur.firstChild.nextSibling && cur.firstChild.nextSibling.type=='textarea')
	divText=cur.firstChild.nextSibling;
      else {
	if (cur.firstChild.type=='textarea')
	  divText=cur.firstChild;
	else return;
      }
    }
    UI.clear_focus();
    if (!UI.readonly)
      UI.scroll_controls();
    // auto-save done
    var s,div;
    //console.log(divText);
    var isdiv=divText.tagName=='DIV',istextarea=divText.type=='textarea',istextinput=divText.tagName =='TEXTINPUT';
    var iscm=!isdiv && !istextarea && !istextinput;
    if (iscm){
      UI.tableaucm=-1;
      s=divText.getValue();
      let f=divText.getTextArea();
      div=f.parentNode;
      UI.cmtotextarea(divText);
      divText=f;
      //iscm=false;
      // console.log('iscm',s);
    }
    else if (istextarea){
      divText.blur();
      s=divText.value;
      divText.innerHTML=s;
      div=divText.parentNode;
    }
    else if (isdiv){
      divText.blur();
      s=divText.innerHTML; // should hide the tex source code instead of replacing?
      s=UI.eval_cmdline1cb(s,s);
      //console.log(s);
      divText.innerHTML=s;
      return;
    }
    else {
      divText.blur();
      s=divText.innerHTML;
      div=divText.parentNode;
    }
    //console.log(divText.classList);
    if (divText.classList.contains('xcas') || divText.classList.contains('cas'))
      UI.ckswitch('xcas',0);
    if (divText.classList.contains('python') || divText.classList.contains('py'))
      UI.ckswitch('python',0);
    if (divText.classList.contains('javascript') || divText.classList.contains('js'))
      UI.ckswitch('javascript',0);
    UI.tableaulocalsave('tableaunoirxcas');
    // remove result of previous evaluation
    //console.log(div);
    var field=div.lastChild;
    if (0 && iscm){
      console.log('iscm',div);
    }
    else {
      for (;field && field!=divText;field=field.previousSibling)
	div.removeChild(field);
    }
    // print
    var divprint=document.createElement("div");
    divprint.style.fontsize=div.style.fontsize;
    divprint.style.width = BoardManager.scrollQuantity()*.9+"px";
    //console.log(divprint.style);
    div.appendChild(divprint);
    // eval
    var s1=s.replace(/<br>/g, '\n');
    s1=s1.replace(/&nbsp;/g,' ');
    if (UI.debug)
      console.log('tableau.js ctrl-enter ',UI.micropy,s1);
    var s2=UI.ckswitch(s1);
    if (s2!=-2)
      return s2;
    // call evaluator
    Module.print_target=divprint;
    try {
      if (UI.micropy==1)
	s2=UI.mp_eval(s1);
      else if (UI.micropy==-1){
 	s2=UI.quickjs(s1);
      }
      else if (UI.micropy==0){
	s2=UI.caseval(s1);
	s2=UI.eval_cmdline1cb(s2,s1);
      }
    } catch (error){
      s2='error: '+error;
    }
    let sp=divprint.innerHTML,dn=0; // remove last end of line
    console.log('s2',s2);
    if (s2.length>=6 && s2[0]=='"' && (s2.substr(1,4)=='<svg' || s2.substr(1,5)=='gr2d(')){
      s2=s2.substr(1,s2.length-2);
    }
    if (s2.length>5 && s2.substr(0,4)=='<svg' && sp.length>9 && sp.substr(sp.length-9,9)=='<br></em>'){
      divprint.innerHTML=sp.substr(0,sp.length-9)+'</em>';
      dn=1;
    }
    //console.log(sp);
    let is_3d= (s2.length>9 && s2.substr(5,4)=='gl3d'),
	is_3d_=(s2.length>=4 && s2.substr(0,4)=='gl3d');
    var gr3ds='';
    if (is_3d || is_3d_){
      var n3d = is_3d?s2.substr(10, s2.length - 16):s2.substr(5,s2.length-5);
      gr3ds='gl3d_'+n3d;
      console.log('gl3d',n3d);
      s2='<div> <canvas id="gl3d_' + n3d + '" onmousedown="UI.canvas_pushed=true;UI.canvas_lastx=event.clientX; UI.canvas_lasty=event.clientY;" onmouseup="UI.canvas_pushed=false;" onmousemove="UI.canvas_mousemove(event,' + n3d + ')" width=' + 400 + ' height=' + 250 + '></canvas></div>';
      //console.log('create3d',s2);
      $id('canvas3d').style.display='none';
    }
    else {
      // translate svg plot
      // console.log(s2);
      if (s2.length>5 && (s2.substr(0,4)=='<svg' || s2.substr(0,5)=='gr2d(' ) )
	s2='<div>'+s2+'</div>';
      // if (s2.length>9 && s2.substr(0,9)=='<div>gl3d') s2='3d not yet supported';
      if (s2.length>9 && s2.substr(0,9)=='<div><svg'){
	//console.log(divText);
	//console.log(divprint);
	//console.log(s1);
	let n=1+UI.count_newline(s1)+UI.count_br(divprint.innerHTML)+dn;
	//console.log('svg n',n);
	n=(4+parseInt(div.style.fontSize))*n+'px"';
	s2='<div><svg style="margin-top:'+n+s2.substr(9,s2.length-9);
      }
      if (s2.length>11 && s2.substr(0,10)=='<div>gr2d('){
	var gr2ds;
	for (;++UI.gr2d_ncanvas;){
	  gr2ds = 'gr2d_' + UI.gr2d_ncanvas;
	  if ($id(gr2ds)==null)
	    break;
	}
	s2 = "<div style='text-align:center'> <canvas id='" + gr2ds + "' width=" + UI.canvas_w + " height=" + UI.canvas_h + "> </canvas><textarea style='display:none'>" + s2.substr(5,s2.length-11) + "</textarea></div>";
	//console.log(s2);
	UI.gr2d_ncanvas++;
      }
    }
    //console.log(s1,' js => ',s2);
    s=divprint.innerHTML;
    s=s+'<!-- -->&nbsp;&nbsp;&nbsp;&nbsp;'+s2;
    //console.log(s);
    divprint.innerHTML=s;
    if (is_3d){
      var el=$id(gr3ds);
      el.addEventListener("touchstart", UI.touch_handler, false);
      el.addEventListener("touchend", UI.touch_handler, false);
      el.addEventListener("touchmove", UI.touch_handler, false);
    }
    Module.print_target=0;
    UI.render_canvas(divprint);
    //divText.id=div.id+'txt';
    if (UI.usemathjax){
      eval("MathJax.typeset();")
      //console.log(divText.id);
      //MathJax.Hub.Process(divText.id);
    }
    // var rect = divText.getBoundingClientRect();
    // console.log(rect.right,rect.top,rect.bottom);
    //var m=MagnetManager.addMagnetText(rect.right+3,rect.bottom+3,s2);
    //m.blur();
    window.getSelection().removeAllRanges();
  } ,
  rewritestring: function (s) {
    var res, i, l, ch;
    l = s.length;
    res = '';
    for (i = 0; i < l; ++i) {
      ch = s[i];
      if (ch == '&') {
        res += "%26";
        continue;
      }
      if (ch == '#') {
        res += "%23";
        continue;
      }
      if (i < l - 2 && ch == '%' && s[i + 1] != '7') {
        res += "___";
        continue;
      }
      res += s[i];
    }
    return res;
  },
  copystringtoclipboard:function(str){
    // Create new element
    var el = document.createElement('textarea');
    // Set value (string to be copied)
    el.value = str;
    // Set non-editable to avoid focus and move outside of view
    el.setAttribute('readonly', '');
    el.style = {position: 'absolute', left: '-9999px'};
    document.body.appendChild(el);
    // Select text inside element
    el.select();
    // Copy text to clipboard
    document.execCommand('copy');
    // Remove temporary element
    document.body.removeChild(el);
  },    
  not_empty:function(txt){
    var s=txt.length;
    //console.log(s);
    for (var i=0;i<s;++i){
      if (txt[i]!=' ' && txt[i]!='\n')
	return true;
    }
    return false;
  },
  is_alphan: function (c) {
    return (c >= 48 && c <= 57) || (c >= 65 && c <= 91) || (c >= 97 && c <= 123) || c == 95;
  },
  python_mode_str: function(i){
    if (i==0) return 'xcas';
    if (i==1) return 'pyth xor';
    if (i==2) return 'pyth **';
    if (i & 4) return 'Python';
    return '?';
  },
  eval_cmdline1cb: function (out, value) {
    var s,graphe=false;
    //console.log('eval_cmdline1cb',out);
    if (out.length > 5 && (out.substr(1, 4) == '<svg' || out.substr(0, 5) == 'gl3d ' || out.substr(0, 5) == 'gr2d(')) {
      //console.log(s);
      s = out;
      out = 'Done_graphic';
      graphe=true;
    }
    else {
      if (out.length > 1 && out[out.length - 1] == ';')
        out = out.substr(0, out.length - 1);
      if (out[0] == '"' || UI.micropy)
        s = 'text ' + out;
      else {
        if (UI.prettyprint && out.length<8192) {
          if (UI.usemathjax)
            s = 'latex(quote(' + out + '))';
          else
            s = 'mathml(quote(' + out + '),1)'; //Module.print(s);
          if (UI.debug)
	    console.log('eval_cmdline1cb ',out,s);
          if (out.length > 10 && out.substr(0, 10) == 'GIAC_ERROR')
            s = '"' + out.substr(11, out.length - 11) + '"';
          else s = UI.caseval_noautosimp(s);
        } else s = out;
      }
    }
    if (s[0]=='"' && s.length>=2) s=s.substr(1,s.length-2);
    if (graphe) s='<div>'+s+'</div>';
    else {
      if (UI.usemathjax) s='\\['+s+'\\]';
    }
    return s;//UI.eval_cmdline1end(value, out, s);
  },
  remove_3d_canvas: function(field){
    var n = field.id;
    if (n && n.length > 5 && n.substr(0, 5) == 'gl3d_') {
      field.parentNode.removeChild(field);
      return;
    }
    var f = field.firstChild;
    for (; f; f = f.nextSibling) {
      UI.remove_3d_canvas(f);
    }
  },
  render_canvas: function (field) {
    // return; // does not work,
    var n = field.id;
    //console.log('render_canvas',n);
    if (n && n.length > 5 && n.substr(0, 5) == 'gr2d_') {
      //console.log(field.nextSibling.value);
      UI.turtle_draw(n, field.nextSibling.value);
    }
    if (n && n.length > 5 && n.substr(0, 5) == 'gl3d_') {
      field.onpointerdown = (e) => { e.stopPropagation(); }
      field.onpointermove = (e) => { e.stopPropagation(); }
      field.onpointerup = (e) => { e.stopPropagation(); }
      var n3d = n.substr(5, n.length - 5);
      console.log(n3d);
      UI.giac_renderer(n3d);
      //Module.canvas=$id('canvas');
      //return;
    }
    var f = field.firstChild;
    for (; f; f = f.nextSibling) {
      UI.render_canvas(f);
    }
  },
  count_newline: function (text) {
    var k = 0, r = 0;
    for (; k < text.length; ++k) {
      if (text.charCodeAt(k) == 10)
        ++r;
    }
    //console.log('count_newline',r);
    return r;
  },
  count_br: function(s){
    //console.log('count_br',s);
    let n=0,pos=0;
    for (;;){
      pos=s.indexOf('<br>');
      if (pos<0 || pos>=s.length)
	return n;
      ++n;
      s=s.substr(pos+4,s.length-pos-4);
    }
  },
  xcascmd: ["!","!=","#","$","%","%/","%/","%{%}","&&","&*","&^","'","()","*","*=","+","+&","+=","+infinity","-","-<","-=","->","-infinity",".*",".+",".-","./",".^","/%","/=",":=","<","<=","=","=<","==","=>",">",">=","?","?","@","@@","ACOSH","ACOT","ACSC","ASEC","ASIN","ASINH","ATAN","ATANH","Airy_Ai","Airy_Bi","Archive","BesselJ","BesselY","Beta","BlockDiagonal","COND","COS","COSH","COT","CSC","CST","Celsius2Fahrenheit","Ci","Circle","ClrDraw","ClrGraph","ClrIO","Col","CopyVar","CyclePic","DIGITS","DOM_COMPLEX","DOM_FLOAT","DOM_FUNC","DOM_IDENT","DOM_INT","DOM_LIST","DOM_RAT","DOM_STRING","DOM_SYMBOLIC","DOM_int","DelFold","DelVar","Det","Dialog","Digits","Dirac","Disp","DispG","DispHome","DrawFunc","DrawInv","DrawParm","DrawPol","DrawSlp","DropDown","DrwCtour","ERROR","EXP","Ei","EndDlog","FALSE","Factor","Fahrenheit2Celsius","False","Fill","GF","Gamma","Gcd","GetFold","Graph","Heaviside","IFTE","Input","InputStr","Int","Inverse","JordanBlock","LN","LQ","LSQ","LU","LambertW","Li","Line","LineHorz","LineTan","LineVert","NORMALD","NewFold","NewPic","Nullspace","Output","Ox_2d_unit_vector","Ox_3d_unit_vector","Oy_2d_unit_vector","Oy_3d_unit_vector","Oz_3d_unit_vector","Pause","Phi","Pi","PopUp","Psi","QR","Quo","REDIM","REPLACE","RandSeed","RclPic","Rem","Request","Resultant","Row","RplcPic","Rref","SCALE","SCALEADD","SCHUR","SIN","SVD","SVL","SWAPCOL","SWAPROW","SetFold","Si","SortA","SortD","StoPic","Store","TAN","TRUE","TeX","Text","Title","True","UTPC","UTPF","UTPN","UTPT","Unarchiv","VARS","VAS","VAS_positive","WAIT","Zeta","[..]","^","_(cm/s)","_(ft/s)","_(m/s)","_(m/s^2)","_(rad/s)","_(rad/s^2)","_(tr/min)","_(tr/s)","_A","_Angstrom","_Bq","_Btu","_Ci","_F","_F_","_Fdy","_G_","_Gal","_Gy","_H","_Hz","_I0_","_J","_K","_Kcal","_MHz","_MW","_MeV","_N","_NA_","_Ohm","_P","_PSun_","_Pa","_R","_REarth_","_RSun_","_R_","_Rankine","_Rinfinity_","_S","_St","_StdP_","_StdT_","_Sv","_T","_V","_Vm_","_W","_Wb","_Wh","_a","_a0_","_acre","_alpha_","_angl_","_arcmin","_arcs","_atm","_au","_b","_bar","_bbl","_bblep","_bu","_buUS","_c3_","_c_","_cal","_cd","_chain","_cm","_cm^2","_cm^3","_ct","_cu","_d","_dB","_deg","_degreeF","_dyn","_eV","_epsilon0_","_epsilon0q_","_epsilonox_","_epsilonsi_","_erg","_f0_","_fath","_fbm","_fc","_fermi","_flam","_fm","_ft","_ft*lb","_ftUS","_ft^2","_ft^3","_g","_g_","_ga","_galC","_galUK","_galUS","_gf","_gmol","_gon","_grad","_grain","_h","_h_","_ha","_hbar_","_hp","_in","_inH20","_inHg","_in^2","_in^3","_j","_kWh","_k_","_kg","_kip","_km","_km^2","_knot","_kph","_kq_","_l","_lam","_lambda0_","_lambdac_","_lb","_lbf","_lbmol","_lbt","_lep","_liqpt","_lm","_lx","_lyr","_m","_mEarth_","_m^2","_m^3","_me_","_mho","_miUS","_miUS^2","_mi^2","_mil","_mile","_mille","_ml","_mm","_mmHg","_mn","_mol","_mp_","_mph","_mpme_","_mu0_","_muB_","_muN_","_oz","_ozUK","_ozfl","_ozt","_pc","_pdl","_ph","_phi_","_pk","_psi","_ptUK","_q_","_qe_","_qepsilon0_","_qme_","_qt","_rad","_rad_","_rd","_rem","_rod","_rpm","_s","_s","_sb","_sd_","_sigma_","_slug","_sr","_st","_syr_","_t","_tbsp","_tec","_tep","_tex","_therm","_ton","_tonUK","_torr","_tr","_tsp","_twopi_","_u","_yd","_yd^2","_yd^3","_yr","_µ","_Âµ","a2q","abcuv","about","abs","abscissa","accumulate_head_tail","acos","acos2asin","acos2atan","acosh","acot","acsc","acyclic","add","add_arc","add_edge","add_vertex","additionally","adjacency_matrix","adjoint_matrix","affix","algsubs","algvar","all_trig_solutions","allpairs_distance","alog10","alors","altitude","and","angle","angle_radian","angleat","angleatraw","animate","animate3d","animation","ans","antiprism_graph","append","apply","approx","approx_mode","arc","arcLen","arccos","arccosh","archive","arclen","arcsin","arcsinh","arctan","arctanh","area","areaat","areaatraw","areaplot","arg","args","array","arrivals","articulation_points","as_function_of","asc","asec","asin","asin2acos","asin2atan","asinh","assert","assign","assign_edge_weights","assume","at","atan","atan2acos","atan2asin","atanh","atrig2ln","augment","auto_correlation","autosimplify","avance","avgRC","axes","back","backquote","backward","baisse_crayon","bandwidth","bar_plot","bartlett_hann_window","barycenter","base","basis","batons","begin","bellman_ford","bernoulli","besselJ","besselY","betad","betad_cdf","betad_icdf","betavariate","bezier","bezout_entiers","biconnected_components","binomial","binomial_cdf","binomial_icdf","bins","bipartite","bipartite_matching","bisection_solver","bisector","bit_depth","bitand","bitor","bitxor","black","blackman_harris_window","blackman_window","bloc","blockmatrix","blue","bohman_window","border","boxwhisker","break","breakpoint","brent_solver","bvpsolve","by","c1oc2","c1op2","cFactor","cSolve","cZeros","cache_tortue","camembert","canonical_form","canonical_labeling","cap","cap_flat_line","cap_round_line","cap_square_line","cartesian_product","cas_setup","case","cat","catch","cauchy","cauchy_cdf","cauchy_icdf","cauchyd","cauchyd_cdf","cauchyd_icdf","cd","cdf","ceil","ceiling","center","center2interval","centered_cube","centered_tetrahedron","cfactor","cfsolve","changebase","channel_data","channels","char","charpoly","chinrem","chisquare","chisquare_cdf","chisquare_icdf","chisquared","chisquared_cdf","chisquared_icdf","chisquaret","choice","cholesky","choosebox","chr","chrem","chromatic_index","chromatic_number","chromatic_polynomial","circle","circumcircle","classes","clear","click","clique_cover","clique_cover_number","clique_number","clique_stats","close","clustering_coefficient","coeff","coeffs","col","colDim","colNorm","colSwap","coldim","collect","colnorm","color","colspace","colswap","comDenom","comb","combine","comment","common_perpendicular","companion","compare","complete_binary_tree","complete_graph","complete_kary_tree","complex","complex_mode","complex_variables","complexroot","concat","cond","cone","confrac","conic","conj","conjugate_equation","conjugate_gradient","connected","connected_components","cont","contains","content","continue","contourplot","contract_edge","convert","convertir","convex","convexhull","convolution","coordinates","copy","correlation","cos","cos2sintan","cosh","cosine_window","cot","cote","count","count_eq","count_inf","count_sup","courbe_parametrique","courbe_polaire","covariance","covariance_correlation","cpartfrac","crationalroot","crayon","createwav","cross","crossP","cross_correlation","cross_point","cross_ratio","crossproduct","csc","csolve","csv2gen","cube","cumSum","cumsum","cumulated_frequencies","curl","current_sheet","curvature","curve","cyan","cycle2perm","cycle_graph","cycleinv","cycles2permu","cyclotomic","cylinder","dash_line","dashdot_line","dashdotdot_line","dayofweek","dayofweek","de","deSolve","debug","debut_enregistrement","default","degree","degree_sequence","del","delcols","delete_arc","delete_edge","delete_vertex","delrows","deltalist","denom","densityplot","departures","derive","deriver","desolve","dessine_tortue","det","det_minor","developper","developper_transcendant","dfc","dfc2f","diag","diff","digraph","dijkstra","dim","directed","discard_edge_attribute","discard_graph_attribute","discard_vertex_attribute","disjoint_union","display","disque","disque_centre","distance","distance2","distanceat","distanceatraw","div","divergence","divide","divis","division_point","divisors","divmod","divpc","dnewton_solver","do","dodecahedron","domain","dot","dotP","dot_paper","dotprod","double","draw_arc","draw_circle","draw_graph","draw_line","draw_pixel","draw_polygon","draw_rectangle","droit","droite_tangente","dsolve","duration","e","e2r","ecart_type","ecart_type_population","ecm_factor","ecris","edge_connectivity","edges","efface","egcd","egv","egvl","eigVc","eigVl","eigenvals","eigenvalues","eigenvectors","eigenvects","element","elif","eliminate","ellipse","else","end","end_for","end_if","end_while","entry","envelope","epaisseur","epaisseur_ligne_1","epaisseur_ligne_2","epaisseur_ligne_3","epaisseur_ligne_4","epaisseur_ligne_5","epaisseur_ligne_6","epaisseur_ligne_7","epaisseur_point_1","epaisseur_point_2","epaisseur_point_3","epaisseur_point_4","epaisseur_point_5","epaisseur_point_6","epaisseur_point_7","epsilon","epsilon2zero","equal","equal2diff","equal2list","equation","equilateral_triangle","erase","erase3d","erf","erfc","error","est_permu","et","euler","euler_gamma","euler_lagrange","eval","eval_level","evala","evalb","evalc","evalf","evalm","even","evolute","exact","exbisector","excircle","execute","exp","exp2list","exp2pow","exp2trig","expand","expexpand","expln","exponential","exponential_cdf","exponential_icdf","exponential_regression","exponential_regression_plot","exponentiald","exponentiald_cdf","exponentiald_icdf","export_graph","expovariate","expr","expression","extend","extract_measure","extrema","ezgcd","f2nd","fMax","fMin","fPart","faces","facteurs_premiers","factor","factor_xn","factorial","factoriser","factoriser_entier","factoriser_sur_C","factors","fadeev","faire","false","falsepos_solver","fclose","fcoeff","fdistrib","feuille","ffaire","ffonction","fft","ffunction","fi","fieldplot","filled","fin_enregistrement","find","findhelp","fisher","fisher_cdf","fisher_icdf","fisherd","fisherd_cdf","fisherd_icdf","fitdistr","flatten","float","float2rational","floor","flow_polynomial","fmod","foldl","foldr","fonction	","fonction_derivee","fopen","for","format","forward","fourier_an","fourier_bn","fourier_cn","fpour","fprint","frac","fracmod","frame_2d","frame_3d","frames","frequencies","frobenius_norm","from","froot","fsi","fsi","fsolve","ftantque","fullparfrac","func","funcplot","function","function_diff","fxnd","gammad","gammad_cdf","gammad_icdf","gammavariate","gauche","gauss","gauss15","gauss_seidel_linsolve","gaussian_window","gaussjord","gaussquad","gbasis","gbasis_max_pairs","gbasis_simult_primes","gcd","gcdex","genpoly","geometric","geometric_cdf","geometric_icdf","getDenom","getKey","getNum","getType","get_edge_attribute","get_edge_weight","get_graph_attribute","get_vertex_attribute","girth","gl_ortho","gl_quaternion","gl_rotation","gl_showaxes","gl_shownames","gl_texture","gl_x","gl_x_axis_color","gl_x_axis_name","gl_x_axis_unit","gl_xtick","gl_y","gl_y_axis_color","gl_y_axis_name","gl_y_axis_unit","gl_ytick","gl_z","gl_z_axis_color","gl_z_axis_name","gl_z_axis_unit","gl_ztick","gnuplot","goto","grad","gramschmidt","graph","graph2tex","graph3d2tex","graph_automorphisms","graph_charpoly","graph_complement","graph_diameter","graph_equal","graph_join","graph_power","graph_rank","graph_spectrum","graph_union","graph_vertices","graphe","graphe3d","graphe_suite","greduce","greedy_color","green","grid_graph","grid_paper","groupermu","hadamard","half_cone","half_line","halftan","halftan_hyp2exp","halt","hamdist","hamming_window","hann_poisson_window","hann_window","harmonic_conjugate","harmonic_division","has","has_arc","has_edge","hasard","head","heading","heapify","heappop","heappush","hermite","hessenberg","hessian","heugcd","hexagon","hidden_name","highlight_edges","highlight_subgraph","highlight_trail","highlight_vertex","highpass","hilbert","histogram","hold","homothety","horner","hybrid_solver","hybridj_solver","hybrids_solver","hybridsj_solver","hyp2exp","hyperbola","hypercube_graph","i","iPart","i[]","iabcuv","ibasis","ibpdv","ibpu","icdf","ichinrem","ichrem","icomp","icontent","icosahedron","id","identifier","identity","idivis","idn","iegcd","if","ifactor","ifactors","ifft","ifte","igamma","igcd","igcdex","ihermite","ilaplace","im","imag","image","implicitdiff","implicitplot","import_graph","in","inString","in_ideal","incidence_matrix","incident_edges","incircle","increasing_power","independence_number","indets","index","induced_subgraph","inequationplot","inf","infinity","infnorm","input","inputform","insert","insmod","int","intDiv","integer","integrate","integrer","inter","interactive_odeplot","interactive_plotode","interp","intersect","interval","interval2center","interval_graph","inv","inverse","inversion","invisible_point","invlaplace","invztrans","iquo","iquorem","iratrecon","irem","isPrime","is_acyclic","is_arborescence","is_biconnected","is_bipartite","is_clique","is_collinear","is_concyclic","is_conjugate","is_connected","is_coplanar","is_cospheric","is_cut_set","is_cycle","is_directed","is_element","is_equilateral","is_eulerian","is_forest","is_graphic_sequence","is_hamiltonian","is_harmonic","is_harmonic_circle_bundle","is_harmonic_line_bundle","is_included","is_inside","is_integer_graph","is_isomorphic","is_isosceles","is_network","is_orthogonal","is_parallel","is_parallelogram","is_permu","is_perpendicular","is_planar","is_prime","is_pseudoprime","is_rectangle","is_regular","is_rhombus","is_square","is_strongly_connected","is_strongly_regular","is_tournament","is_tree","is_triconnected","is_two_edge_connected","is_vertex_colorable","is_weighted","ismith","isobarycenter","isom","isomorphic_copy","isopolygon","isosceles_triangle","isprime","ithprime","jacobi_equation","jacobi_linsolve","jacobi_symbol","jordan","jusqu_a","jusqua","jusque","kde","keep_algext","keep_pivot","ker","kernel","kernel_density","kill","kneser_graph","kolmogorovd","kolmogorovt","kovacicsols","l1norm","l2norm","label","labels","lagrange","laguerre","laplace","laplacian","laplacian_matrix","latex","lcf_graph","lcm","lcoeff","ldegree","left","left_rectangle","legend","legendre","legendre_symbol","len","length","leve_crayon","lgcd","lhs","ligne_chapeau_carre","ligne_chapeau_plat","ligne_chapeau_rond","ligne_polygonale","ligne_polygonale_pointee","ligne_tiret","ligne_tiret_point","ligne_tiret_pointpoint","ligne_trait_plein","limit","limite","lin","line","line_graph","line_inter","line_paper","line_segments","line_width_1","line_width_2","line_width_3","line_width_4","line_width_5","line_width_6","line_width_7","linear_interpolate","linear_regression","linear_regression_plot","lineariser","lineariser_trigo","linfnorm","linsolve","linspace","lis","lis_phrase","list","list2exp","list2mat","list_edge_attributes","list_graph_attributes","list_vertex_attributes","listplot","lll","ln","lname","lncollect","lnexpand","local","locus","log","log10","logarithmic_regression","logarithmic_regression_plot","logb","logistic_regression","logistic_regression_plot","lower","lowest_common_ancestor","lowpass","lp_assume","lp_bestprojection","lp_binary","lp_binaryvariables","lp_breadthfirst","lp_depthfirst","lp_depthlimit","lp_firstfractional","lp_gaptolerance","lp_hybrid","lp_initialpoint","lp_integer","lp_integertolerance","lp_integervariables","lp_interiorpoint","lp_iterationlimit","lp_lastfractional","lp_maxcuts","lp_maximize","lp_method","lp_mostfractional","lp_nodelimit","lp_nodeselect","lp_nonnegative","lp_nonnegint","lp_pseudocost","lp_simplex","lp_timelimit","lp_variables","lp_varselect","lp_verbose","lpsolve","lsmod","lsq","lu","lvar","mRow","mRowAdd","magenta","make_directed","make_weighted","makelist","makemat","makesuite","makevector","map","maple2mupad","maple2xcas","maple_ifactors","maple_mode","markov","mat2list","mathml","matpow","matrix","matrix_norm","max","maxflow","maximal_independent_set","maximize","maximum_clique","maximum_degree","maximum_independent_set","maximum_matching","maxnorm","mean","median","median_line","member","mgf","mid","middle_point","midpoint","min","minimal_edge_coloring","minimal_spanning_tree","minimal_vertex_coloring","minimax","minimize","minimum_degree","minus","mkisom","mksa","mod","modgcd","mods","montre_tortue","moustache","moyal","moyenne","mul","mult_c_conjugate","mult_conjugate","multinomial","multiplier_conjugue","multiplier_conjugue_complexe","multiply","mupad2maple","mupad2xcas","mycielski","nCr","nDeriv","nInt","nPr","nSolve","ncols","negbinomial","negbinomial_cdf","negbinomial_icdf","neighbors","network_transitivity","newList","newMat","newton","newton_solver","newtonj_solver","nextperm","nextprime","nlpsolve","nodisp","nom_cache","non","non_recursive_normal","nop","nops","norm","normal","normal_cdf","normal_icdf","normald","normald_cdf","normald_icdf","normalize","normalt","normalvariate","not","nprimes","nrows","nuage_points","nullspace","number_of_edges","number_of_spanning_trees","number_of_triangles","number_of_vertices","numer","octahedron","od","odd","odd_girth","odd_graph","odeplot","odesolve","of","op","open","open_polygon","option","or","ord","order","order_size","ordinate","orthocenter","orthogonal","osculating_circle","otherwise","ou","output","p1oc2","p1op2","pa2b2","pade","parabola","parallel","parallelepiped","parallelogram","parameq","parameter","paramplot","parfrac","pari","part","partfrac","parzen_window","pas","pas_de_cote","path_graph","pcar","pcar_hessenberg","pcoef","pcoeff","pencolor","pendown","penup","perimeter","perimeterat","perimeteratraw","periodic","perm","perminv","permu2cycles","permu2mat","permuorder","permute_vertices","perpen_bisector","perpendicular","petersen_graph","peval","pi","piecewise","pivot","pixoff","pixon","planar","plane","plane_dual","playsnd","plex","plot","plot3d","plotarea","plotcdf","plotcontour","plotdensity","plotfield","plotfunc","plotimplicit","plotinequation","plotlist","plotode","plotparam","plotpolar","plotproba","plotseq","plotspectrum","plotwav","plus_point","pmin","point","point2d","point3d","point_carre","point_croix","point_etoile","point_invisible","point_losange","point_milieu","point_plus","point_point","point_triangle","point_width_1","point_width_2","point_width_3","point_width_4","point_width_5","point_width_6","point_width_7","poisson","poisson_cdf","poisson_icdf","poisson_window","polar","polar_coordinates","polar_point","polarplot","pole","poly2symb","polyEval","polygon","polygone_rempli","polygonplot","polygonscatterplot","polyhedron","polynom","polynomial_regression","polynomial_regression_plot","position","poslbdLMQ","posubLMQ","potential","pour","pow","pow2exp","power_regression","power_regression_plot","powermod","powerpc","powexpand","powmod","prepend","preval","prevperm","prevprime","primpart","print","printf","prism","prism_graph","proc","product","program","projection","proot","propFrac","propfrac","psrgcd","ptayl","purge","pwd","pyramid","python_compat","q2a","qr","quadrant1","quadrant2","quadrant3","quadrant4","quadric","quadrilateral","quantile","quartile1","quartile3","quartiles","quest","quo","quorem","quote","r2e","radical_axis","radius","ramene","rand","randMat","randNorm","randPoly","randbetad","randbinomial","randchisquare","randexp","randfisher","randgammad","randgeometric","randint","randmarkov","randmatrix","randmultinomial","randnorm","random","random_bipartite_graph","random_digraph","random_graph","random_network","random_planar_graph","random_regular_graph","random_sequence_graph","random_tournament","random_tree","random_variable","randperm","randpoisson","randpoly","randseed","randstudent","randvar","randvector","randweibulld","range","rank","ranm","ranv","rassembler_trigo","rat_jordan","rational","rationalroot","ratnormal","rcl","rdiv","re","read","readrgb","readwav","real","realroot","reciprocation","rectangle","rectangle_droit","rectangle_gauche","rectangle_plein","rectangular_coordinates","recule","red","redim","reduced_conic","reduced_quadric","ref","reflection","regroup","relabel_vertices","reliability_polynomial","rem","remain","remove","reorder","repeat","repete","repeter","replace","resample","residue","resoudre","resoudre_dans_C","resoudre_systeme_lineaire","restart","resultant","return","reverse","reverse_graph","reverse_rsolve","revert","revlex","revlist","rgb","rhombus","rhombus_point","rhs","riemann_window","right","right_rectangle","right_triangle","risch","rm_a_z","rm_all_vars","rmbreakpoint","rmmod","rmwatch","romberg","rombergm","rombergt","rond","root","rootof","roots","rotate","rotation","round","row","rowAdd","rowDim","rowNorm","rowSwap","rowdim","rownorm","rowspace","rowswap","rref","rsolve","same","sample","samplerate","sans_factoriser","saute","sauve","save_history","scalarProduct","scalar_product","scale","scaleadd","scatterplot","schur","sec","secant_solver","segment","seidel_spectrum","seidel_switch","select","semi_augment","seq","seqplot","seqsolve","sequence_graph","series","set[]","set_edge_attribute","set_edge_weight","set_graph_attribute","set_pixel","set_vertex_attribute","set_vertex_positions","shift","shift_phase","shortest_path","show_pixels","shuffle","si","sierpinski_graph","sign","signature","signe","similarity","simp2","simplex_reduce","simplifier","simplify","simpson","simult","sin","sin2costan","sincos","single_inter","sinh","sinon","size","sizes","slope","slopeat","slopeatraw","smith","smith","smod","snedecor","snedecor_cdf","snedecor_icdf","snedecord","snedecord_cdf","snedecord_icdf","solid_line","solve","somme","sommet","sort","sorta","sortd","sorted","soundsec","spanning_tree","sphere","spline","split","spring","sq","sqrfree","sqrt","square","square_point","srand","sst","sst_in","st_ordering","stack","star_graph","star_point","start","stdDev","stddev","stddevp","steffenson_solver","step","stereo2mono","sto","str","string","string","strongly_connected_components","student","student_cdf","student_icdf","studentd","studentt","sturm","sturmab","sturmseq","style","subMat","subdivide_edges","subgraph","subs","subsop","subst","substituer","subtype","sum","sum_riemann","suppress","surd","svd","swapcol","swaprow","switch","switch_axes","sylvester","symb2poly","symbol","syst2mat","tCollect","tExpand","table","tablefunc","tableseq","tabvar","tail","tan","tan2cossin2","tan2sincos","tan2sincos2","tangent","tangente","tanh","tantque","taux_accroissement","taylor","tchebyshev1","tchebyshev2","tcoeff","tcollect","tdeg","tensor_product","test","tetrahedron","texpand","textinput","then","thickness","thiele","threshold","throw","time","title","titre","tlin","to","tonnetz","topologic_sort","topological_sort","torus_grid_graph","tourne_droite","tourne_gauche","tpsolve","trace","trail","trail2edges","trames","tran","transitive_closure","translation","transpose","trapeze","trapezoid","traveling_salesman","tree","tree_height","triangle","triangle_paper","triangle_plein","triangle_point","triangle_window","trig2exp","trigcos","trigexpand","triginterp","trigsimplify","trigsin","trigtan","trn","true","trunc","truncate","try","tsimplify","tuer","tukey_window","tutte_polynomial","two_edge_connected_components","type","ufactor","ugamma","unapply","unarchive","underlying_graph","unfactored","uniform","uniform_cdf","uniform_icdf","uniformd","uniformd_cdf","uniformd_icdf","union","unitV","unquote","until","upper","user_operator","usimplify","valuation","vandermonde","var","variables_are_files","variance","vector","vector","vers","version","vertex_connectivity","vertex_degree","vertex_distance","vertex_in_degree","vertex_out_degree","vertices","vertices_abc","vertices_abca","vpotential","watch","web_graph","weibull","weibull_cdf","weibull_icdf","weibulld","weibulld_cdf","weibulld_icdf","weibullvariate","weight_matrix","weighted","weights","welch_window","wheel_graph","when","while","white","widget_size","wilcoxonp","wilcoxons","wilcoxont","with_sqrt","write","writergb","writewav","wz_certificate","xcas_mode","xor","xyztrange","yellow","zeros","zip","ztrans","{}","|","||"],
  dicho_find: function (tableau, s) {
    var l = tableau.length, debut = 0, fin = l, milieu;
    if (l == 0) return false;
    if (s < tableau[0] || s > tableau[l - 1]) return false;
    // s>=tableau[debut] and s<=tableau[fin-1]
    for (; debut < fin - 1;) {
      milieu = Math.floor((debut + fin) / 2);
      // console.log(debut,fin,milieu,tableau[milieu])
      if (s >= tableau[milieu]) debut = milieu; else fin = milieu;
    }
    // console.log(s,tableau[debut]);
    if (s == tableau[debut]) return true;
    return false;
  },
  unique_completion: function (tableau, s) {
    var l = tableau.length, debut = 0, fin = l, milieu;
    if (l == 0 || s > tableau[l - 1]) return "";
    if (s < tableau[0]) return s == tableau[0].substr(s.length);
    // s>=tableau[debut] and s<=tableau[fin-1]
    for (; debut < fin - 1;) {
      milieu = Math.floor((debut + fin) / 2);
      //console.log(debut,fin,milieu,tableau[milieu])
      if (s >= tableau[milieu]) debut = milieu; else fin = milieu;
    }
    // console.log(tableau[debut+1],tableau[debut+1].length);
    // console.log(s,tableau[debut],tableau[debut+1],tableau[debut+2]);
    // s>=tableau[debut]
    if (s == tableau[debut]) {
      //if (debut + 1 < l && s == tableau[debut + 1].substr(0, s.length)) return "";
      return s;
    }
    // now s>tableau[debut] and s is supposed to be shorter, hence s!=begin of tableau[debut]
    if (debut + 1 < l && s != tableau[debut + 1].substr(0, s.length)) return "";
    if (debut + 2 < l && s == tableau[debut + 2].substr(0, s.length)) return "";
    return tableau[debut + 1];
  },
  renderhelp: function (text, localized_cmd) {
    var s = text.substr(1, text.length - 2);
    var pos0 = s.search("</b>");
    var found = (s.substr(pos0 + 7, 20) != "Best match has score");
    var lh = s.substr(3, pos0 - 3);
    var unlocalized_cmd = lh;
    //console.log(localized_cmd,unlocalized_cmd);
    if (found) {
      var tmp = eval('longhelp.' + lh);
      if (tmp == undefined) lh = '';
      else {
        var lang = tmp.substr(7, 2);
        var prefix = UI.docprefix;
        prefix = prefix.substr(0, prefix.length - 13) + lang + '/cascmd_' + lang + '/';
        //console.log(prefix,tmp);
        lh = ' (<a href="' + prefix + tmp + '" target="_blank">' + (UI.langue == -1 ? '+ de d&eacute;tails' : 'more details') + '</a>)';
      }
    }
    //Module.print(lh);
    var sorig = s;
    pos1 = s.search("<br>");
    if (pos1 < 0) return sorig;
    var explication = s.substr(0, pos1);
    s = s.substr(pos1 + 4, s.length - pos1 - 4);
    pos1 = s.search("<br>");
    if (pos1 < 0) return sorig;
    var syntaxe = s.substr(0, pos1);
    s = s.substr(pos1 + 4, s.length - pos1 - 4);
    pos1 = s.search("<br>");
    if (pos1 < 0) return sorig;
    var voiraussi = s.substr(0, pos1);
    var examples = s.substr(pos1 + 4, s.length - pos1 - 4);
    if (found)
      s = explication + lh + '<br><tt>' + syntaxe + '</tt><br>' + (UI.langue == -1 ? 'Voir aussi: ' : 'See also ');
    else
      s = explication.substr(0, pos0);
    while (true) {
      pos1 = voiraussi.search(',');
      if (pos1 < 0) break;
      var cmd = UI.clean(voiraussi.substr(0, pos1), false);
      if (found)
        s += '<button class="bouton" onmousedown="event.preventDefault()" onclick="UI.addhelp(\'?\',\'' + cmd + '\')">' + cmd + '</button>';
      else
        s += '<button class="bouton" onmousedown="event.preventDefault()" onclick="var help=$id(\'helptxt\');help.value=\'' + cmd + '\'; UI.insert_focused(\'' + UI.clean(cmd, true) + '(\')">' + cmd + '</button>';
      voiraussi = voiraussi.substr(pos1 + 1, voiraussi.length - pos1 - 1);
    }
    if (found)
      s += '<button class="bouton" onmousedown="event.preventDefault()" onclick="UI.addhelp(\'?\',\'' + voiraussi + '\')">' + voiraussi + '</button>';
    else
      s += '<button class="bouton" onmousedown="event.preventDefault()" onclick="UI.insert_focused(\'' + voiraussi + '(\')">' + voiraussi + '</button>';
    pos1 = examples.search(';');
    if (pos1 >= 0) {
      s += '<br>' + (UI.langue == -1 ? 'Exemples: ' : 'Examples ');
      while (true) {
        pos1 = examples.search(';');
        if (pos1 < 0) break;
        var cmd = examples.substr(0, pos1);
        //console.log('avant',cmd);
        cmd = cmd.replace(unlocalized_cmd, localized_cmd);
        //console.log('apres',cmd);
        cmd = UI.clean(cmd, false);
        s += '<button class="bouton" onmousedown="event.preventDefault()" onclick="UI.insert_focused(\'' + UI.clean(cmd, true) + '\')">' + cmd + '</button>';
        examples = examples.substr(pos1 + 1, examples.length - pos1 - 1);
      }
    }
    s += '<button class="bouton" onmousedown="event.preventDefault()" onclick="UI.insert_focused(\'' + UI.clean(examples, true) + '\')">' + examples + '</button>';
    //console.log(s);
    return s;
  },
  addhelp: function (prefixe, text) {
    $id('online_help').style.display='block';
    $id('helptxt').value = text;
    var input = prefixe + text;
    var out = UI.caseval(input);
    var add=UI.renderhelp(out,text);
    var helpoutput = $id('helpoutput');
    helpoutput.innerHTML = add;
    //console.log(helpoutput);
  },
  completion: function (cm,showhint=true) {
    var s, k;
    if (cm.type == 'textarea') {
      k = cm.selectionStart;
      s = cm.value;
    }
    else {
      var pos = cm.getCursor();
      k = pos.ch;
      s = cm.getLine(pos.line);
    }
    var kstart = k;
    // skip at cursor
    for (; k > 0; k--) {
      var c = s.charCodeAt(k);
      if (UI.is_alphan(c))
        break;
    }
    var kend = k;
    for (; k >= 0; k--) {
      var c = s.charCodeAt(k);
      if (!UI.is_alphan(c))
        break;
    }
    for (; k < kend; k++) {
      var c = s.charCodeAt(k + 1);
      if (c > 64) break;
    }
    //Module.print(s); Module.print(k); Module.print(kend);
    if (s.length < 2) {
      $id('online_help').style.display='block';
      //UI.insert(UI.focused, '?');
      return;
    }
    s = s.substr(k + 1, kend - k);
    var sc = UI.unique_completion(UI.xcascmd, s);
    //console.log('completion',s,sc,showhint);
    //Module.print(s); Module.print(sc);
    if (cm.type == 'textarea') {
      cm.selectionStart = k + 1;
      cm.selectionEnd = kstart;
      if (sc != "") {
        UI.insert(cm, sc);
        s = sc;
        cm.selectionStart = k + 1;
        cm.selectionEnd = k + 1 + sc.length;
      }
    } else {
      //console.log(sc.length);
      if (sc != "") {
        //console.log(s,k);
        cm.setSelection({line: pos.line, ch: k + 1}, {line: pos.line, ch: pos.ch});
        UI.insert(cm, sc);
        //cm.setSelection({line:pos.line,ch:k+1},{line:pos.line,ch:k+1+sc.length});
        cm.setSelection({line: pos.line, ch: k + 1 + sc.length}, {line: pos.line, ch: k + 1 + sc.length});
        s = sc;
      }
      else {
        //console.log(s,k);
	if (showhint){
          cm.setCursor({line: pos.line, ch: k + 1 + s.length});
          cm.showHint();
	  return;
	}
      }
    }
    UI.addhelp('?', s);
  },
  isie: false,
  scrollatend: function (field) {
    if (!UI.isie)
      field.scrollTop = field.scrollHeight;
  },
  getsel: function (field) {
    var startPos = field.selectionStart;
    var endPos = field.selectionEnd;
    var selectedText = field.value.substring(startPos, endPos);
    return selectedText;
  },
  indent_or_complete: function (field) {
    if (field.type != 'textarea') {
      if (field.lineCount() == 1)
        UI.completion(field);
      else
        field.execCommand('indentAuto');
    }
  },
  indentline: function (field) {
    if (field.type != 'textarea') {
      field.execCommand('indentAuto');
    }
  },
  moveCaret: function (field, charCount) {
    if (field.type != "textarea" && field.type != "text") {
      var pos = field.getCursor();
      pos.ch = pos.ch + charCount;
      field.setCursor(pos);
      field.refresh();
      // UI.show_curseur();
      return;
    }
    var pos = field.selectionStart;
    pos = pos + charCount;
    if (pos < 0) pos = 0;
    if (pos > field.value.length) pos = field.value.length;
    field.setSelectionRange(pos, pos);
  },
  show_curseur: function () {
    document.getElementsByClassName("CodeMirror-cursors")[0].style.visibility = "visible";
    var cursors = document.getElementsByClassName("CodeMirror-cursor");
    for (var i = 0; i < cursors.length; i++) {
      cursors[i].style.visibility = "visible";
    }
  },
  moveCaretUpDown: function (field, Count) {
    if (field.type != "textarea" && field.type != "text") {
      var pos = field.getCursor();
      pos.line = pos.line + Count;
      field.setCursor(pos);
      field.refresh();
      //UI.show_curseur();
      return;
    }
    if (Count < -1) {
      var i;
      for (i = 0; i > Count; i--)
        UI.moveCaretUpDown(field, -1);
      return;
    }
    if (Count > 1) {
      var i;
      for (i = 0; i < Count; i++)
        UI.moveCaretUpDown(field, 1);
      return;
    }
    var pos = field.selectionStart;
    var s = field.value;
    var cur = pos, shift = pos + 1, pos1;
    cur--;
    if (cur >= s.length) cur--;
    for (; cur >= 0; cur--) {
      if (s.charCodeAt(cur) == 10) {
        shift = pos - cur;
        break;
      }
    }
    if (Count == -1) {
      if (cur < 0) return;
      pos1 = cur;
      cur--;
      for (; cur >= 0; cur--) {
        if (s.charCodeAt(cur) == 10) break;
      }
      //console.log(cur,shift);
      pos = cur + shift;
      if (pos > pos1) pos = pos1;
    }
    if (Count == 1) {
      cur = pos;
      for (; cur < s.length; cur++) {
        if (s.charCodeAt(cur) == 10) break;
      }
      pos = cur + shift;
      if (pos >= s.length) return;
      pos1 = pos;
      for (; pos1 > cur; pos1--) {
        if (s.charCodeAt(pos1) == 10) pos = pos1;
      }
    }
    if (pos < 0) pos = 0;
    if (pos > field.value.length) pos = field.value.length;
    field.setSelectionRange(pos, pos);
  },
  backspace: function (field) {
    //if (UI.focusaftereval) field.focus();
    if (field.type != "textarea" && field.type != "text") {
      var start = field.getCursor('from');
      var end = field.getCursor('to');
      if (end.line != start.line || end.ch != start.ch)
        field.replaceSelection('');
      else {
        var c = start.ch;
        var l = start.line;
        if (start.ch == 0 && start.line == 0) return;
        if (c > 0) {
          var s = field.getRange({line: l, ch: 0}, end), i;
          for (i = 0; i < s.length; i++) {
            if (s.charAt(i) != ' ') break;
          }
          //console.log(i,s.length,c);
          if (i == s.length && c >= 2) {
            var l1 = l - 1;
            for (; l1 >= 0; --l1) {
              s = field.getLine(l1);
              for (i = 0; i < s.length && i < c; i++) {
                if (s[i] != ' ') break;
              }
              if (i != s.length && i < c) break;
            }
            if (l1 >= 0) c = i; else c -= 2;
          }
          else c--;
          field.replaceRange('', {line: l, ch: c}, end);
        }
        else {
          l--;
          var s = field.getRange({line: l, ch: c}, end);
          field.replaceRange('', {line: l, ch: s.length - 1}, end);
        }
      }
      var t = field.getTextArea();
      t.value = field.getValue();
    } else {
      var pos = field.selectionStart;
      var pos2 = field.selectionEnd;
      var s = field.value;
      if (pos < pos2) {
        field.value = s.substring(0, pos) + s.substring(pos2, s.length);
        if (pos < 0) pos = 0;
        if (pos > field.value.length) pos = field.value.length;
        field.setSelectionRange(pos, pos);
        UI.resizetextarea(field);
        return;
      }
      if (pos > 0) {
        field.value = s.substring(0, pos - 1) + s.substring(pos, s.length);
        pos--;
        if (pos < 0) pos = 0;
        if (pos > field.value.length) pos = field.value.length;
        field.setSelectionRange(pos, pos);
        UI.resizetextarea(field);
      }
    }
  },
  resizetextarea: function (field) {
    if (field.type != 'textarea') return;
    var s = field.value;
    var C = field.cols;
    //console.log(C);
    var N = 0, i, j = 0, n = s.length, c;
    for (i = 0; i < n; i++, j++) {
      c = s.charCodeAt(i);
      if (c == 10 || j == C) {
        N++;
        j = 0;
      }
    }
    field.rows = N + 1;
  },
  changefontsize: function (field, size,maxw=0,maxh=0) {
    let f=field.getWrapperElement();
    f.style["font-size"] = size + "px";
    if (maxw!=0) // does not scroll
      f.style["max-width"]=maxw+"px";
    if (maxh!=0)
      f.style["max-height"]=maxh+"px";
    field.refresh();
  },
  // logo turtle and pixel library display
  color_list: ['black',
    'red',
    'green',
    'yellow',
    'blue',
    'magenta',
    'cyan',
    'white',
    'silver',
    'gray',
    'maroon',
    'purple',
    'fuchsia',
    'lime',
    'olive',
    'navy',
    'teal',
    'aqua',
    'antiquewhite',
    'aquamarine',
    'azure',
    'beige',
    'bisque',
    'blanchedalmond',
    'blueviolet',
    'brown',
    'burlywood',
    'cadetblue',
    'chartreuse',
    'chocolate',
    'coral',
    'cornflowerblue',
    'cornsilk',
    'crimson',
    'cyan',
    'darkblue',
    'darkcyan',
    'darkgoldenrod',
    'darkgray',
    'darkgreen',
    'darkgrey',
    'darkkhaki',
    'darkmagenta',
    'darkolivegreen',
    'darkorange',
    'darkorchid',
    'darkred',
    'darksalmon',
    'darkseagreen',
    'darkslateblue',
    'darkslategray',
    'darkslategrey',
    'darkturquoise',
    'darkviolet',
    'deeppink',
    'deepskyblue',
    'dimgray',
    'dimgrey',
    'dodgerblue',
    'firebrick',
    'floralwhite',
    'forestgreen',
    'gainsboro',
    'ghostwhite',
    'gold',
    'goldenrod',
    'greenyellow',
    'grey',
    'honeydew',
    'hotpink',
    'indianred',
    'indigo',
    'ivory',
    'khaki',
    'lavender',
    'lavenderblush',
    'lawngreen',
    'lemonchiffon',
    'lightblue',
    'lightcoral',
    'lightcyan',
    'lightgoldenrodyellow',
    'lightgray',
    'lightgreen',
    'lightgrey',
    'lightpink',
    'lightsalmon',
    'lightseagreen',
    'lightskyblue',
    'lightslategray',
    'lightslategrey',
    'lightsteelblue',
    'lightyellow',
    'limegreen',
    'linen',
    'mediumaquamarine',
    'mediumblue',
    'mediumorchid',
    'mediumpurple',
    'mediumseagreen',
    'mediumslateblue',
    'mediumspringgreen',
    'mediumturquoise',
    'mediumvioletred',
    'midnightblue',
    'mintcream',
    'mistyrose',
    'moccasin',
    'navajowhite',
    'oldlace',
    'olivedrab',
    'orangered',
    'orchid',
    'palegoldenrod',
    'palegreen',
    'paleturquoise',
    'palevioletred',
    'papayawhip',
    'peachpuff',
    'peru',
    'pink',
    'plum',
    'powderblue',
    'rosybrown',
    'royalblue',
    'saddlebrown',
    'salmon',
    'sandybrown',
    'seagreen',
    'seashell',
    'sienna',
    'skyblue',
    'slateblue',
    'slategray',
    'slategrey',
    'snow',
    'springgreen',
    'steelblue',
    'tan',
    'thistle',
    'tomato',
    'turquoise',
    'violet',
    'wheat',
    'whitesmoke',
    'yellowgreen'],
  arc_en_ciel: function (k) {
    var r, g, b;
    k += 21;
    k %= 126;
    if (k < 0)
      k += 126;
    if (k < 21) {
      r = 251;
      g = 0;
      b = 12 * k;
    }
    if (k >= 21 && k < 42) {
      r = 251 - (12 * (k - 21));
      g = 0;
      b = 251;
    }
    if (k >= 42 && k < 63) {
      r = 0;
      g = (k - 42) * 12;
      b = 251;
    }
    if (k >= 63 && k < 84) {
      r = 0;
      g = 251;
      b = 251 - (k - 63) * 12;
    }
    if (k >= 84 && k < 105) {
      r = (k - 84) * 12;
      g = 251;
      b = 0;
    }
    if (k >= 105 && k < 126) {
      r = 251;
      g = 251 - (k - 105) * 12;
      b = 0;
    }
    return 'rgb(' + r + ',' + g + ',' + b + ')';
  },
  turtle_color: function (c) {
    if (c >= 0x100) {
      if (c < 0x17e)
        return UI.arc_en_ciel(c);
      //console.log('rgb('+Math.floor(c/(256*256))+','+(Math.floor(c/256) % 256)+','+(c%256)+')');
      var r=8*((c>>11) & 0x1f);
      var g=4*((c>>5) & 0x3f);
      var b=8*(c & 0x1f);
      return 'rgb(' + r + ',' + g + ',' + b + ')';
      // return 'rgb(' + Math.floor(c / (256 * 256)) + ',' + (Math.floor(c / 256) % 256) + ',' + (c % 256) + ')';
    }
    return UI.color_list[c];
  },
  pixon_draw: function (id, s) {
    var v = eval(s);
    if (!Array.isArray(v)) return;
    //console.log(v[0], v.length);
    var canvas = $id(id);
    var l = v.length, w = 0, h = 0;
    if (l < 2) return;
    var scale = v[0];
    for (var k = 1; k < l; k++) {
      var cur = v[k];
      var x = cur[0], y = cur[1];
      if (cur.length==3 && typeof cur[2]!="number"){
	x+=100;
	y+=16;
      }
      if (cur.length==4) {
        var tmp = cur[3];
	if (typeof tmp=="number"){
          if (tmp > 0) y += tmp; else x -= tmp;
	} else {
	  x+=100;
	  y+=16;
	}
      }
      //console.log(cur,x,y);
      if (x > w) w = x;
      if (y > h) h = y;
    }
    w = (w + 1) * scale;
    h = (h + 1) * scale;
    canvas.width = w;
    canvas.height = h;
    //console.log(h,w);
    if (canvas.getContext) {
      var ctx = canvas.getContext('2d');
      for (var k = 1; k < l; k++) {
        var cur = v[k], cl;
        //console.log(cur);
        if (!Array.isArray(cur) || (cl = cur.length) < 2) continue;
        // cur[0]=x, cur[1]=y, cur[2]=color, cur[3]=w if +, h if -
        var x = cur[0] * scale;
        var y = cur[1] * scale;
	if (cl>2 && typeof cur[2]=="string"){
	  console.log(cur[2]);
	  ctx.font = '16px serif';
	  ctx.fillStyle = 'black';
	  ctx.fillText(cur[2],x,y+16,100);
	  continue;
	}
        ctx.fillStyle = (cl > 2) ? UI.turtle_color(cur[2]) : 'black';
        if (cl < 4) {
          ctx.fillRect(x, y, scale, scale);
          continue;
        }
	if (typeof cur[3]=="string"){
	  ctx.font = '16px serif';
	  ctx.fillText(cur[3],x,y+16,100);
	  continue;
	}
        var h = cur[3] * scale, w = scale;
        if (h < 0) {
          w = -h;
          h = scale;
        }
        ctx.fillRect(x, y, w, h);
      }
    }
  },
  turtle_dx: 0, // shift frame
  turtle_dy: 0,
  turtle_z: 1,  // zoom factor
  turtle_maillage: 1,
  turtle_draw: function (id, s) {
    if (s.length < 7) return;
    s = s.substr(5, s.length - 6);
    //console.log(s);
    if (s.length > 7 && s.substr(s, 6) == "pixon(") {
      UI.pixon_draw(id, s.substr(6, s.length - 7));
      return;
    }
    if (s.length < 6 || s.substr(s, 5) != "logo(")
      return;
    s = s.substr(5, s.length - 6);
    //console.log(s);
    var v = eval(s);
    if (!Array.isArray(v)) return;
    //console.log(v[0]);
    var canvas = $id(id);
    var w = canvas.width, h = canvas.height;
    if (canvas.getContext) {
      var ctx = canvas.getContext('2d');
      var turtlezoom = UI.turtle_z, turtlex = UI.turtle_dx, turtley = UI.turtle_dy;
      // background
      ctx.fillStyle = 'white';
      ctx.fillRect(0,0,w,h);
      // maillage
      if (UI.turtle_maillage & 3) {
        ctx.fillStyle = 'black';
        var xdecal = Math.floor(turtlex / 10.0) * 10;
        var ydecal = Math.floor(turtley / 10.0) * 10;
        if ((UI.turtle_maillage & 0x3) == 1) {
          for (var i = xdecal; i < w + xdecal; i += 10) {
            for (var j = ydecal; j < h + ydecal; j += 10) {
              var X = Math.floor((i - turtlex) * turtlezoom + .5);
              var Y = Math.floor((j - turtley) * turtlezoom + .5);
              // console.log(X,Y);
              ctx.fillRect(X, h - Y, 1, 1);
            }
          }
        } else {
          var dj = Math.sqrt(3.0) * 10, i0 = xdecal;
          for (var j = ydecal; j < h + ydecal; j += dj) {
            var J = Math.floor(h - (j - turtley) * turtlezoom);
            for (var i = i0; i < w + xdecal; i += 10) {
              ctx.fillRect(Math.floor((i - turtlex) * turtlezoom + .5), J, 1, 1);
            }
            i0 += dj;
            while (i0 >= 10)
              i0 -= 10;
          }
        }
      }
      var l = v.length, i;
      // montre la position et le cap (v[l-1])
      var prec = v[l - 1];
      ctx.font = '16px serif';
      ctx.fillStyle = 'yellow';
      ctx.fillRect(w - 40, 0, 40, 50);
      ctx.fillStyle = 'black';
      ctx.fillText('x:' + prec[0], w - 40, 15);
      ctx.fillText('y:' + prec[1], w - 40, 31);
      ctx.fillText('t:' + prec[2], w - 40, 49);
      // v[i]=[x(0),y(1),cap(2),status(3),r(4),chaine(5)],
      // couleur=status >> 11
      // longueur_tortue= (status>>3)&0xff
      // direct=status&4 (vrai si angle dans le sens trigo)
      // visible=status&2
      // crayon baisse=status&1
      // si r>0 arc/disque rayon=r & 0x1ff, theta1=(r >> 9) & 0x1ff, theta2=(r >> 18) & 0x1ff
      //        rempli=(r>>27)&0x1
      // si r<0 ligne polygonale extremite v[i] origine v[i+r] (r<0)
      for (k = 1; k < l; k++) {
        prec = v[k - 1];
        var cur = v[k];
        var preccouleur = prec[3] >> 11; // -> FIXME colors
        var curcouleur = prec[3] >> 11; // -> FIXME colors
        if (cur[5].length) {
          ctx.font = cur[4] + 'px serif';
          ctx.strokeStyle = ctx.fillStyle = UI.turtle_color(curcouleur);
          ctx.fillText(cur[5], turtlezoom * (cur[0] - turtlex), h - turtlezoom * (cur[1] - turtley));
          continue;
        }
        var radius = cur[4], precradius = prec[4];
        var x1 = Math.floor(turtlezoom * (prec[0] - turtlex) + .5),
            y1 = Math.floor(turtlezoom * (prec[1] - turtley) + .5),
            x2 = Math.floor(turtlezoom * (cur[0] - turtlex) + .5),
            y2 = Math.floor(turtlezoom * (cur[1] - turtley) + .5);
        if (radius > 0) {
          var r = radius & 0x1ff, theta1, theta2, rempli, x, y, R, angle;
          theta1 = prec[2]+ ((radius >> 9) & 0x1ff);
          theta2 = prec[2] + ((radius >> 18) & 0x1ff);
          rempli = (radius >> 27) & 1;
	  var seg = (radius >> 28) & 1;
          R = Math.floor(turtlezoom * r + .5);
          angle1 = Math.PI / 180 * (theta1 - 90);
          angle2 = Math.PI / 180 * (theta2 - 90);
          x = Math.floor(turtlezoom * (cur[0] - turtlex - r * Math.cos(angle2)) + .5);
          y = Math.floor(turtlezoom * (cur[1] - turtley - r * Math.sin(angle2)) + .5);
          ctx.beginPath();
	  if (seg)
            ctx.moveTo(x2, h - y2);
	  else {
            ctx.moveTo(x, h - y);
            ctx.lineTo(x2, h - y2);
	  }
	  //console.log(x,y,x1,y1,angle1,angle2);
          ctx.arc(x, h - y, R, -angle2,-angle1);
          ctx.closePath();
          ctx.strokeStyle = ctx.fillStyle = UI.turtle_color(curcouleur);
          if (rempli)
            ctx.fill();
          else
            ctx.stroke();
          continue;
        }
        if (prec[3] & 1) {
          ctx.strokeStyle = ctx.fillStyle = UI.turtle_color(preccouleur);
          ctx.beginPath();
          ctx.moveTo(x1, h - y1);
          ctx.lineTo(x2, h - y2);
          ctx.closePath();
          ctx.stroke();
        }
        if (radius < -1 && k + radius >= 0) {
          ctx.strokeStyle = ctx.fillStyle = UI.turtle_color(curcouleur);
          ctx.beginPath();
          var x0 = Math.floor(turtlezoom * (cur[0] - turtlex) + .5), y0 = Math.floor(turtlezoom * (cur[1] - turtley) + .5);
          //console.log('begin',x0,y0);
          ctx.moveTo(x0, h - y0);
          for (var i = -1; i >= radius; i--) {
            prec = v[k + i];
            var x = Math.floor(turtlezoom * (prec[0] - turtlex) + .5);
            var y = Math.floor(turtlezoom * (prec[1] - turtley) + .5);
            //console.log(i,x,y);
            ctx.lineTo(x, h - y);
          }
          //console.log('end',x0,y0);
          //ctx.lineTo(x0,h-y0);
          ctx.closePath();
          ctx.fill(); // automatically close path
        }
      }
      var cur = v[l - 1];
      if (cur[3] & 2) {
        // dessin de la tortue
        var x = Math.floor(turtlezoom * (cur[0] - turtlex) + .5);
        var y = Math.floor(turtlezoom * (cur[1] - turtley) + .5);
        var cost = Math.cos(cur[2] * Math.PI / 180);
        var sint = Math.sin(cur[2] * Math.PI / 180);
        var turtle_length = (cur[3] >> 3) & 0xff;
        var Dx = Math.floor(turtlezoom * turtle_length * cost / 2 + .5);
        var Dy = Math.floor(turtlezoom * turtle_length * sint / 2 + .5);
        //console.log('tortue',cur,w,h,turtlezoom,x,y,Dx,Dy);
        ctx.strokeStyle = ctx.fillStyle = UI.turtle_color(curcouleur);
        ctx.beginPath();
        ctx.moveTo(x + Dy, h - (y - Dx));
        ctx.lineTo(x - Dy, h - (y + Dx));
        ctx.closePath();
        ctx.stroke();
        if (!(cur[3] & 1))
          ctx.strokeStyle = ctx.fillStyle = UI.turtle_color(curcouleur + 1);
        ctx.beginPath();
        ctx.moveTo(x + Dy, h - (y - Dx));
        ctx.lineTo(x + 3 * Dx, h - (y + 3 * Dy));
        ctx.closePath();
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(x - Dy, h - (y + Dx));
        ctx.lineTo(x + 3 * Dx, h - (y + 3 * Dy));
        ctx.closePath();
        ctx.stroke();
      }
    }
  },
  p2inp1p3:function(p1x,p1y,p2x,p2y,p3x,p3y,eps){
    // check if p2(p2x,p2y) is inside [p1,p3] up to precision eps
    let dx2=p2x-p1x,dy2=p2y-p1y,dx3=p3x-p1x,dy3=p3y-p1y;
    let d2=dx2*dx2+dy2*dy2,d3=dx3*dx3+dy3*dy3;
    if (d2==0 || (p2x==p3x && p2y==p3y))
      return true;
    // distance test p3 must be farther than p2 from p1
    if (d3<d2) return false;
    // test by angle cosine
    let c=(dx2*dx3+dy2*dy3);
    if (c<=0) return false;
    c=(c*c)/(d2*d3);
    // must be near 1
    if (1-c<=eps*eps/d2){
      //if (c<1) console.log('p2inp1p3',p1x,p1y,p2x,p2y,p3x,p3y);
      return true;
    }
    return false;
  },
  encode:function(points,pressures){
    let l=points.length;
    if (l!=pressures.length){
      console.log('error UI.encode points.length!=pressures.length');
      return '';
    }
    if (l<2)
      return ''; // point
    let lastp=-1,lastx=points[0].x,lasty=points[0].y,res='';
    for (let pos=1;pos<l;++pos){
      if (pressures[pos]!=lastp){ 
	res += pressures[pos]; // always if pos==1
	lastp=pressures[pos];
      }
      let dx=points[pos].x-lastx;
      lastx=points[pos].x;
      let dy=points[pos].y-lasty;
      lasty=points[pos].y;
      res += UI.relative(dx,1)+UI.relative(dy,0);
      //console.log(dx,dy,res);
    }
    //console.log(res);
    return res;
  },
  drawline:function(context, x1, y1, x2, y2, pressure = 1.0, color = user.getCurrentColor()) {
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
  },
  drawdot:function(context,x, y, color) {
    context.beginPath();
    context.fillStyle = color;
    context.lineWidth = 2.5;
    context.arc(x, y, 2, 0, 2 * Math.PI);
    context.fill();
    context.closePath();
  },
  draw_encoded:function(s,ctx=0,DX=0,DY=0){
    if (ctx==0)
      ctx=document.getElementById("canvas").getContext("2d");
    //console.log(s);
    // polyline drawing encoded by a string formatted by
    // color,startX,startY,relative strokes (see UI.relative for relative positions encoding)
    let lxypc=UI.decode_polygon(s,DX,DY);
    let lx=lxypc[0],ly=lxypc[1],lp=lxypc[2],color=lxypc[3];
    let lastX=lx[0],lastY=ly[0],curX,curY,pressure=0.5;
    if (lx.length==1){
      UI.drawdot(ctx,lastX,lastY,color);
      return;
    }
    for (let i=1;i<lx.length;++i){
      curX=lx[i]; curY=ly[i];
      UI.drawline(ctx,lastX,lastY,curX,curY,lp[i],color);
      lastX=curX; lastY=curY;
    }
    /*
    let l=s.length,pos=s.search(','),lastX=0,lastY=0;
    let color="white",pressure=0.5;
    if (s[0]>'9'){
      color=s.substr(0,pos);
      s=s.substr(pos+1,l-pos-1);
      pos=s.search(',');
    }
    lastX=eval(s.substr(0,pos))+DX;
    s=s.substr(pos+1,l-pos-1);
    pos=s.search(',');
    lastY=eval(s.substr(0,pos))+DY;
    // console.log(color,lastX,lastY); return;
    l=s.length; s=s.substr(pos+1,l-pos-1);
    // now s has format [digit]majuscules minuscules
    l=s.length;
    if (l==0){
      //console.log(color,lastX,lastY,pressure);
      drawDot(lastX,lastY,color); return;
    }
    for (pos=0;pos<l;){
      let n=s.charCodeAt(pos);
      if (n>=48 && n<=57){
	pressure=(n-48)/10+0.05;
	//console.log('pressure',pressure);
	++pos;
      }
      let dx=0,dy=0;
      for (;pos<l;++pos){
	n=s.charCodeAt(pos);
	if (n>96)
	  break;
	if (n<65){
	  console.log('invalid char');
	  return;
	}
	dx *= 13;
	dx += (n-77); // 65+12 where 'A'==65
      }
      // console.log(dx);
      for (;pos<l;++pos){
	n=s.charCodeAt(pos);
	if (n<97)
	  break;
	if (n>122){
	  console.log('invalid char');
	  return;
	}
	dy *= 13;
	dy += (n-109); // 97+12 where 'A'==65
      }
      // console.log(dy);
      let curX=lastX+dx;
      let curY=lastY+dy;
      //console.log(dx,dy);
      drawLine(ctx,lastX,lastY,curX,curY,pressure,color);
      lastX=curX; lastY=curY;
    }
    */
  },
  relative:function(dx,majuscule){
    dx=Math.round(dx);
    // encode dx with A-Z: M=0, A-L -12..-1, N..Y 1..12
    let tabM=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y'];
    let tabm=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y'];
    if (dx>=-12 && dx<=12)
      return majuscule?tabM[dx+12]:tabm[dx+12];
    // console.log(dx);
    let res='';
    while (dx!=0){
      let q=Math.round(dx/13);
      dx -= q*13;
      res = (majuscule?tabM[dx+12]:tabm[dx+12])+res;
      dx = q;
    }
    //console.log(res);
    return res;
  },
  order_magnets:function(m1,m2){
    let x1=parseInt(m1.style.left),y1=parseInt(m1.style.top),
	x2=parseInt(m2.style.left),y2=parseInt(m2.style.top);
    if (x1==x2)
      return Math.sign(y1-y2);
    return Math.sign(x1-x2);
  },
  reorder_magnets:function(){
    let magnets = document.getElementsByClassName("magnet");
    let tab=[];
    for (let i = 0; i < magnets.length; i++) {
      magnets[i].style.left=Math.floor(parseInt(magnets[i].style.left)/200)*200+'px';
      tab.push(magnets[i]);
    }
    tab.sort(UI.order_magnets);
    let listm=$id("magnets");
    while (listm.hasChildNodes())
      listm.removeChild(listm.firstChild);
    // console.log(tab);
    for (let i=0;i<tab.length;++i){
      listm.appendChild(tab[i]);
    }
    //console.log(listm.children.length);
    UI.link();
  },
  exec_magnets:function(){
    let magnets = document.getElementsByClassName("magnet");
    for (let i = 0; i < magnets.length; i++) {
      UI.tableaueval(magnets[i]);
    }
    UI.link();
  },
  display_trash_msg:true,
  ck_display_trash_msg:function(){
    UI.link();
    if (UI.display_trash_msg){
      window.alert(UI.langue==-1?'Les magnets peuvent etre restaures depuis Menu':'Magnets may be restored from Menu');
      UI.display_trash_msg=false;
    }
  },
  restore_trash:function(){
    let magnets = document.getElementsByClassName("magnet");
    for (let i = 0; i < magnets.length; i++) {
      var cur=magnets[i];
      //console.log('restore',i,cur.style.display);
      if (cur.style.visibility=='hidden')
	cur.style.visibility='visible';
    }    
  },
  empty_trash:function(){
    let magnets = document.getElementsByClassName("magnet");
    for (let i = 0; i < magnets.length; i++) {
      var cur=magnets[i];
      if (cur.style.visibility=='hidden')
	cur.remove();
    }    
    UI.link();
  },
  magnets_to_trash:function(){
    let magnets = document.getElementsByClassName("magnet");
    for (let i = 0; i < magnets.length; i++) {
      magnets[i].style.visibility='hidden';
    }
    UI.ck_display_trash_msg();
  },
  init_xcas:function(){
    $id('readonly_no').style.fontWeight="900";
    if (UI.detectmob()){
      $id("pc_controls").style.display='none';//UI.usecm=false;
    }
    else
      $id("mobile_controls").style.display='none';
    //history.pushState(null, null, window.location.pathname);
    //window.addEventListener('popstate',UI.popstate, false);
    window.onbeforeunload = function (e) {
      var dialogText = UI.langue==-1?"Etes-vous sur?":"Are you sure?";
      return dialogText;
    };
    var ck = UI.readCookie('xcas_magnet_button');
    if (ck){
      UI.ckswitch(ck,0);
    }
    ck = UI.readCookie('xcas_magnet_allow_exec');
    //console.log(ck);
    if (ck=='false') UI.allow(false); else UI.allow(true);
    ck = UI.readCookie('xcas_angle_radian');
    if (ck=='0') UI.radian(false); else UI.radian(true);
    ck = UI.readCookie('xcas_forum_url');
    if (ck){
      UI.forum_url=ck;
      $id('forum_url').value=ck;
    }
    ck = UI.readCookie('xcas_warn_link');
    UI.set_warnlink(ck!='0');
    UI.set_assistants_color();
    svgwidth=8.;
    $id("controls").scrollIntoView();
  },
  scroll_magnets:function(dx){
    // console.log('scroll magnets',dx);
    let magnets = document.getElementsByClassName("magnet");
    //console.log('right',BoardManager.scrollQuantity());
    for (let i = 0; i < magnets.length; i++) {
      // console.log('scroll', magnets[i].style.left,-dx);
      var tmp=parseInt(magnets[i].style.left);
      tmp -= dx;
      magnets[i].style.left = tmp+"px";
    }
  },
  tableau_width:function(){
    // min width required for magnets and delineations
    let w=0,dx=container.scrollLeft;
    let magnets = document.getElementsByClassName("magnet");
    for (let i = 0; i < magnets.length; i++) {
      let tmp=parseInt(magnets[i].style.left)+dx;
      if (magnets[i].style.width)
	tmp += parseInt(magnets[i].style.width);
      else
	tmp += 400;
      // console.log('tableau_width magnet',i,tmp,dx);
      if (tmp>w)
	w=tmp;
    }
    let s=BoardManager.cancelStack.stack,W,H,x,y;
    [W,H,x,y]=UI.find_whxy(s);
    if (x+W>w)
      w=x+W;
    console.log('tableau_width',w);
    return w;
  },
  canvas2png: function(canvasElement, fileName) {
    // save canvas to png
    var MIME_TYPE = "image/png";
    var imgURL = canvasElement.toDataURL(MIME_TYPE);
    var dlLink = document.createElement('a');
    dlLink.download = fileName;
    dlLink.href = imgURL;
    dlLink.dataset.downloadurl = [MIME_TYPE, dlLink.download, dlLink.href].join(':');

    document.body.appendChild(dlLink);
    dlLink.click();
    document.body.removeChild(dlLink);
  },
  prepare_png_export:function(xleft=0){
    let x = container.scrollLeft;
    UI.scroll_magnets(xleft-x);
    container.scrollTo({top:0,left:xleft});
    if (1){
      const nodeContent = document.getElementById("content");
      const nodeBoard = document.getElementById("board");
      let w=canvas.width; //
      w=UI.tableau_width();
      w=Math.ceil(w/800)*800;//console.log(w);
      let h=canvas.height; // $id('canvasBackground').height
      h=1200;
      nodeContent.style.width = "" + w;
      nodeContent.style.height = "" + h;
      nodeBoard.style.width = "" + w;
      nodeBoard.style.height = "" + h;
      nodeContent.style.transition = "";
      nodeContent.style.transform = "scale(1)";
    }
    return x;
  },
  png_export:function(filename,field_id){ // field_id not supported yet
    var hidden=controls.hidden;
    controls.hidden=true;
    Menu.hide();
    let nodeContent = $id("content"),nodeBoard = $id("board"),x=0,w=1200,dx=1200;
    if (field_id)
      nodeContent=$id(field_id);
    else {
      x=UI.prepare_png_export();
      w=UI.tableau_width();
    }
    if (1) html2canvas(nodeContent).then(canvas => {
      UI.canvas2png(canvas,filename);
      controls.hidden=hidden;
    });
    nodeBoard.style.height=nodeBoard.style.width=nodeContent.style.width=nodeContent.style.heigth=nodeContent.style.transition = "";
    UI.scroll_magnets(x);
    container.scrollTo({top:0,left:x});
    controls.hidden=hidden;
  },
  invert_canvas:function(canvas){
    let hidden=controls.hidden;
    controls.hidden=true;
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
  },
  pdf_export:function(filename,invert=0,field_id=0){ // field_id not supported yet
    var hidden=controls.hidden;
    controls.hidden=true;
    Menu.hide();
    let nodeContent = $id("content"),x=0,nodeBoard = $id("board");
    if (field_id)
      nodeContent=$id(field_id);
    else 
      x=UI.prepare_png_export();
    let W=UI.tableau_width(); // W*=1.25;
    // html2canvas(nodeContent, {  windowWidth: W    });
    nodeContent.style.width=""+$id('canvas').width;
    nodeBoard.style.width=""+$id('canvas').width;
    html2canvas(nodeContent,
		{
		  windowWidth: W,
		  windowHeight: nodeContent.scrollHeight
		}
	       ).then(canvas => {
      const doc = new jspdf.jsPDF('l');
      let firstpage = true;
      const h = 1200;//window.innerHeight;
      const w = 1600;//window.innerWidth;
      W=Math.ceil(W/w)*w;
      for (let x = 0; x < W; x += w) {
        if (!firstpage)
          doc.addPage();
        const canvasPage = document.createElement("canvas");
        //console.log(canvas.height)
        canvasPage.height = h;
        canvasPage.width = w;
        const scaleCorrection = canvas.height / h; //because the html2canvas canvas is of height 1090 instead of 1000. Why?
        canvasPage.getContext("2d").drawImage(canvas, x*scaleCorrection, 0, w*scaleCorrection, h*scaleCorrection, 0, 0, w, h);
        if (invert)
	  UI.invert_canvas(canvasPage);
        const pw = doc.internal.pageSize.getWidth();
        const ph = doc.internal.pageSize.getHeight();
	//console.log('print',scaleCorrection,w,h,pw,ph,w/h,pw/ph);
        if (w/h > pw/ph)
          doc.addImage(canvasPage, 0, 0, pw, pw*h/w);
        else
          doc.addImage(canvasPage, 0, 0, ph*w/h, ph);
        firstpage = false;
      }
      doc.save(filename);
      controls.hidden=hidden;
    });
    nodeBoard.style.width=nodeContent.style.width=nodeContent.style.heigth=nodeContent.style.transition = "";
    UI.scroll_magnets(-x);
    container.scrollTo({top:0,left:x});
  },
}; // closing UI={
