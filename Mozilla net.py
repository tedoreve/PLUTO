import urllib
import re
from bs4 import BeautifulSoup as bs
import codecs
import gzip

#======================读取网页=================================================
def getPage(url):
    headers = ('User-Agent','Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11')
    opener = urllib.request.build_opener()
    opener.addheaders = [headers]
    page = opener.open(url)
    return page.read()
#======================无格式输出===============================================
#print(html)
#======================有格式输出===============================================
#h=getPage('http://www.drao.nrc.ca/')
#h = bs(h, 'html.parser')
#print(h.prettify())
#=======================下载图片================================================
def getImg(html):
    reg = r'src="(.+?\.jpg)" pic_ext'
    imgre = re.compile(reg)
    html=html.decode('utf-8')
    imglist = re.findall(imgre,html)
    x = 0
    for imgurl in imglist:
        urllib.request.urlretrieve(imgurl,'%s.jpg' % x)
        x+=1
    return imglist
#===========================下载pdf============================================
def getPdf(html):
    reg = r'href="(.+?\.pdf)" '
    imgre = re.compile(reg)
    html=html.decode('utf-8')
    imglist = re.findall(imgre,html)
    x = 0
    for imgurl in imglist:
        urllib.request.urlretrieve(imgurl,'%s.pdf' % x)
        x+=1
    return imglist
#========================返回某一网页中包含特定字符的链接========================
def getLink(html):
    reg = r'HREF="(snrs.G.+?\.html)"'
    imgre = re.compile(reg)
    html=html.decode('utf-8')
    imglist = re.findall(imgre,html)
    for imgurl in imglist:
        h=getPage('http://www.mrao.cam.ac.uk/surveys/snrs/'+imgurl)
        h=h.decode('utf-8')
        if ('Distance' in h):
            print('http://www.mrao.cam.ac.uk/surveys/snrs/'+imgurl)
    return imglist
#=============================搜索B站视频名=====================================
def getBili(html):
    reg = r'<title>([\s\S]*)</title>'
    imgre = re.compile(reg)
    html=gzip.decompress(html).decode('utf-8')
    imglist = re.findall(imgre,html)
    return imglist
#==============================================================================
#========================主程序================================================
#h=getPage("http://bt.neu6.edu.cn/")
#result=getLink(h)
n=100
result=[]
for i in range(n):
    try:
        h=getPage("http://www.bilibili.com/video/av1"+str(i))
        result=result+getBili(h)
        if i%2==1:
            print(int((i+1)/2)*'#',int((i+1)/n*100),'%')
        if i%2==0:
            print(int(i/2)*'#',int((i+1)/n*100),'%')
    except:
        if i%2==1:
            print(int((i+1)/2)*'#',int((i+1)/n*100),'%')
        if i%2==0:
            print(int(i/2)*'#',int((i+1)/n*100),'%')
#========================保存网页(不完善)=======================================
#h=h.decode('utf-8')
#file = codecs.open('py.txt', 'w','utf-8')
#file.write(h)
#file.close()
#===========================测试链接============================================
#domain=[
#'baidu.com',\
#'google.ac',\
#'google.ad',\
#'google.ae',\
#'google.com.af',\
#'google.com.ag',\
#'google.com.ai',\
#'google.al',\
#'google.am',\
#'google.co.ao',\
#'google.com.ar',\
#'google.as',\
#'google.at',\
#'google.com.au',\
#'google.az',\
#'google.ba',\
#'google.com.bd',\
#'google.be',\
#'google.bf',\
#'google.bg',\
#'google.com.bh',\
#'google.bi',\
#'google.bj',\
#'google.com.bn',\
#'google.com.bo',\
#'google.com.br',\
#'google.bs',\
#'google.bt',\
#'google.co.bw',\
#'google.by',\
#'google.com.bz',\
#'google.ca',\
#'google.com.kh',\
#'google.cc',\
#'google.cd',\
#'google.cf',\
#'google.cat',\
#'google.cg',\
#'google.ch',\
#'google.ci',\
#'google.co.ck',\
#'google.cl',\
#'google.cm',\
#'google.cn',\
#'google.com.co',\
#'google.co.cr',\
#'google.com.cu',\
#'google.cv',\
#'google.com.cy',\
#'google.cz',\
#'google.de',\
#'google.dj',\
#'google.dk',\
#'google.dm',\
#'google.com.do',\
#'google.dz',\
#'google.com.ec',\
#'google.ee',\
#'google.com.eg',\
#'google.es',\
#'google.com.et',\
#'google.fi',\
#'google.com.fj',\
#'google.fm',\
#'google.fr',\
#'google.ga',\
#'google.ge',\
#'google.gf',\
#'google.gg',\
#'google.com.gh',\
#'google.com.gi',\
#'google.gl',\
#'google.gm',\
#'google.gp',\
#'google.gr',\
#'google.com.gt',\
#'google.gy',\
#'google.com.hk',\
#'google.hn',\
#'google.hr',\
#'google.ht',\
#'google.hu',\
#'google.co.id',\
#'google.iq',\
#'google.ie',\
#'google.co.il',\
#'google.im',\
#'google.co.in',\
#'google.io',\
#'google.is',\
#'google.it',\
#'google.je',\
#'google.com.jm',\
#'google.jo',\
#'google.co.jp',\
#'google.co.ke',\
#'google.ki',\
#'google.kg',\
#'google.co.kr',\
#'google.com.kw',\
#'google.kz',\
#'google.la',\
#'google.com.lb',\
#'google.com.lc',\
#'google.li',\
#'google.lk',\
#'google.co.ls',\
#'google.lt',\
#'google.lu',\
#'google.lv',\
#'google.com.ly',\
#'google.co.ma',\
#'google.md',\
#'google.me',\
#'google.mg',\
#'google.mk',\
#'google.ml',\
#'google.com.mm',\
#'google.mn',\
#'google.ms',\
#'google.com.mt',\
#'google.mu',\
#'google.mv',\
#'google.mw',\
#'google.com.mx',\
#'google.com.my',\
#'google.co.mz',\
#'google.com.na',\
#'google.ne',\
#'google.com.nf',\
#'google.com.ng',\
#'google.com.ni',\
#'google.nl',\
#'google.no',\
#'google.com.np',\
#'google.nr',\
#'google.nu',\
#'google.co.nz',\
#'google.com.om',\
#'google.com.pk',\
#'google.com.pa',\
#'google.com.pe',\
#'google.com.ph',\
#'google.pl',\
#'google.com.pg',\
#'google.pn',\
#'google.com.pr',\
#'google.ps',\
#'google.pt',\
#'google.com.py',\
#'google.com.qa',\
#'google.ro',\
#'google.rs',\
#'google.ru',\
#'google.rw',\
#'google.com.sa',\
#'google.com.sb',\
#'google.sc',\
#'google.se',\
#'google.com.sg',\
#'google.sh',\
#'google.si',\
#'google.sk',\
#'google.com.sl',\
#'google.sn',\
#'google.sm',\
#'google.so',\
#'google.st',\
#'google.sr',\
#'google.com.sv',\
#'google.td',\
#'google.tg',\
#'google.co.th',\
#'google.com.tj',\
#'google.tk',\
#'google.tl',\
#'google.tm',\
#'google.to',\
#'google.tn',\
#'google.com.tr',\
#'google.tt',\
#'google.com.tw',\
#'google.co.tz',\
#'google.com.ua',\
#'google.co.ug',\
#'google.co.uk',\
#'google.com',\
#'google.com.uy',\
#'google.co.uz',\
#'google.com.vc',\
#'google.co.ve',\
#'google.vg',\
#'google.co.vi',\
#'google.com.vn',\
#'google.vu',\
#'google.ws',\
#'google.co.za',\
#'google.co.zm',\
#'google.co.zw',\
#]
#for i in domain:
#    try:
#        h=getPage('https://www.'+i)
#        print(i)
#    except:
#        print(0)
#==============================================================================
