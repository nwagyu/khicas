SUBDIRS = 	en es fr el local pt zh de
DIST_SUBDIRS = 	en es fr el local pt zh de 
XCASHTML = xcasfr.html xcasen.html *.css *.js giacwasm.wasm logo.png giac.tex giacfr.tex test.tex testfr.tex test.html testfr.html algolycee.html graphtheory-user_manual.pdf *.sty simple.html *tableau* magnets html2canvas.js logo.png w3.css jspdf.umd.min.js khicasnw.html epsilon.dfu apps.tar recovery nws*.html nws.js numworks.js shakeys khi*.bin bootloader.bin
EXTRA_DIST = aide_cas README index.html $(XCASHTML) dxcas/* pari/* castex.tex tortue.tex 50-numworks-calculator.rules  # private nwscert.cc -> shakeys

dist_pkgdata_DATA = aide_cas
dist_doc_DATA = README index.html


install-exec-hook: 
	$(install_sh) -d $(DESTDIR)$(prefix)/share/giac/doc/dxcas
	$(install_sh) -d $(DESTDIR)$(prefix)/share/giac/doc/pari
	for dd in $(XCASHTML); do \
	  $(INSTALL_DATA) $(top_srcdir)/doc/$$dd $(DESTDIR)$(prefix)/share/giac/doc/ ; \
	done
	$(INSTALL_DATA) $(top_srcdir)/doc/dxcas/* $(DESTDIR)$(prefix)/share/giac/doc/dxcas
	$(INSTALL_DATA) $(top_srcdir)/doc/pari/* $(DESTDIR)$(prefix)/share/giac/doc/pari
	$(INSTALL_DATA) $(top_srcdir)/doc/index.html $(DESTDIR)$(prefix)/share/giac/doc/
	$(INSTALL_DATA) $(top_srcdir)/doc/README $(DESTDIR)$(prefix)/share/giac/doc/
	ln -sf ../aide_cas $(DESTDIR)$(prefix)/share/giac/doc/aide_cas
#	cp -f epsilon.dfu apps.tar recovery $(DESTDIR)$(prefix)/share/giac/doc/
	cp -f 50-numworks-calculator.rules $(DESTDIR)$(prefix)/share/giac/doc/
#	cp -f shakeys $(DESTDIR)$(prefix)/share/giac/doc/
#	ln -sf $(DESTDIR)$(prefix)/share/giac/doc/en/cascmd_en $(DESTDIR)$(prefix)/share/giac/doc/local/cascmd_local
#	ln -sf $(DESTDIR)$(prefix)/share/giac/doc/en/casinter $(DESTDIR)$(prefix)/share/giac/doc/local/casinter
#	ln -sf $(DESTDIR)$(prefix)/share/giac/doc/en/keywords $(DESTDIR)$(prefix)/share/giac/doc/local/keywords
#	ln -sf $(DESTDIR)$(prefix)/share/giac/doc/en/xcasex $(DESTDIR)$(prefix)/share/giac/doc/local/xcasex
#	ln -sf $(DESTDIR)$(prefix)/share/giac/doc/en/xcasmenu $(DESTDIR)$(prefix)/share/giac/doc/local/xcasmenu

