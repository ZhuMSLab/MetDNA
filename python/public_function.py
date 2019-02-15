import os
import shutil
import re

# URL = 'http://localhost:8024/media'
# URL = 'http://metdna.zhulab.cn/media'  # formal cloud url
# URL = 'http://192.168.201.33/media'  # formal web url
URL = 'http://192.168.201.33:8024/media'  # test web url
MEDIA_ROOT = '/mnt/data/metdna-upload'  # local path in docker

THREADS = 5

CE = {"10": "10", "20": "20",
      "30": "30", u"35\u00b115": "35,15", "40": "40",
      "50": "50", "NCE15-30-45": "30"}

SPECIES = {"Homo sapiens (human)":  "hsa", "Mus musculus (mouse)": "mmu",
           "Rattus norvegicus (rat)": "rat", "Bos taurus (cow)": "bta",
           "Drosophila melanogaster (fruit fly)": "dme",
           "Gallus gallus (chicken)": "gga", "Danio rerio (zebrafish)": "dre",
           "Caenorhabditis elegans (nematode)": "cel",
           "Saccharomyces cerevisiae (yeast)": "sce",
           "Arabidopsis thaliana (thale cress)": "ath",
           "Schistosoma mansoni": "smm",
           "Plasmodum falciparum 3D7 (Malaria)": "pfa",
           "Trypanosoma brucei": "tbr",
           "Escherichia coli K-12 MG1655": "eco",
           "Pseudomonas putida KT2440": "ppu", "Synechococcus elongatus": "syf"}

STAT_METHOD = {'Student t-test': 't', 'Wilcox test': 'wilcox'}

ADJUST_P = {'No': 'FALSE', 'Yes': 'TRUE'}

NAME_MAPPING = {'adjustP': 'correct',
                'caseGroup': 'case.group',
                'ce': 'ce',
                'controlGroup': 'control.group',
                'lcType': 'column',
                'pCutoff': 'p.cutoff',
                'polarity': 'polarity',
                'species': 'species',
                'statMethod': 'uni.test',
                'path': 'path',
                'instrument': 'instrument'}

INSTT = {
    'Sciex TripleTOF': 'SciexTripleTOF',
    'Agilent QTOF': 'AgilentQTOF',
    'Other QTOF': 'OtherQTOF',
    'Thermo Orbitrap (HCD)': 'ThermoOrbitrap'
}

a_dic = {'polarity': u'positive', 'control.group': u'g1',
         'column': u'hilic', 'ce': '10', 'species': 'hsa',
         'uni.test': 't', 'p.cutoff': 0.05, 'threads': 10,
         'path': u'/tmp/metdna-upload\\22cd-7466-4355-ac40\\a75381a1',
         'case.group': u'g2', 'correct': u'No'}

RESULT_FILE_NAME = 'results.tar.gz'

with open('./db_passwd.secret') as f:
    db_para = [i.strip() for i in f.read().split('\n')]

_ = {'host': db_para[0], 'user': 'metdna', 'port': int(db_para[1]),
     'passwd': db_para[2], 'db': 'metdna'}  # log in information

Q = {
    'q_pro_queue':  # project queue
        """
        SELECT * FROM `metDNACore_projectqueue` a 
        WHERE a.`status`='waiting' ORDER BY a.`submit_time`;
        """,
    'q_files':  # upload files related to a particular project
        """
        SELECT * FROM `metDNACore_uploadfile` a 
        WHERE a.`project_id`={id:d};
        """,
    'q_user':
        """
        SELECT * FROM `metDNACore_customuser` a WHERE a.`id`={id:d};
        """,
    'q_pro':  # project
        """
        SELECT * FROM `metDNACore_project` a WHERE a.`id`={id:d};
        """,
    'query2':
        """
        SELECT * FROM `metDNACore_projectqueue` a WHERE a.`project_id`=285;
        """,
    'update_before_analysis':  # update 'start_time' in project queue
        """
        UPDATE `metDNACore_projectqueue` SET `status`='{status:s}', `start_time`=NOW() 
        WHERE `project_id`={project_id:d};
        """,
    'update_after_analysis':  # update 'end_time' and 'status' in project queue
        """
        UPDATE `metDNACore_projectqueue` SET `status`='{status:s}', `end_time`=NOW() 
        WHERE `project_id`={project_id:d};
        """
}


BASE_MSG = """\
<html>
  <head></head>
  <body style="padding: 20px 300px 20px 20px">
    <p>Dear MetDNA user,</p>
    <p>Thank you for using MetDNA!</p>
    <p>Your project <b>{project_name:s}</b> is finished and the analysis results can be downloaded within 7 days via the link:</p>
    <p>{url:s}</p>
    <br>
    <p>The detailed description of the analysis results can be found <a href="{demo_link:s}">here</a>.</p>
    <br>
    <p>If you have any questions, please see the <a href="{faq_link:s}">FAQs of MetDNA</a>, or send email to us (metdna@sioc.ac.cn).</p>
    <br>
    <p>--------------------------------------</p>
    <div>
      <p>Laboratory for Mass Spectrometry and Metabolomics (ZHU LAB)</p>
      <p>Interdisciplinary Research Center on Biology and Chemistry (IRCBC)</p>
      <p>Shanghai Institute of Organic Chemistry (SIOC)</p>
      <p>Chinese Academy of Sciences (CAS)</p>
      <p>26 Qiuyue Road, Pudong, Shanghai, China 201210</p>
      <p>Website: <a href="http://www.zhulab.cn">www.zhulab.cn</a></p>
    </div>
  </body>
</html>
"""


BASE_MSG_ERROR = """\
<html>
  <head></head>
  <body style="padding: 20px 300px 20px 20px">
    <p>Dear MetDNA user,</p>
    <p>Thank you for using MetDNA!</p>
    <p>We are very sorry to send you this email because some errors have occurred 
       during the analysis of your project <b>{project_name:s}.</b></p>
    <p>We will check this error later.</p>
    <br>
    <p>If you have any questions, please see the <a href="{faq_link:s}">FAQs of MetDNA</a>, or send email to us (metdna@sioc.ac.cn).</p>
    <br>
    <p>--------------------------------------</p>
    <div>
      <p>Laboratory for Mass Spectrometry and Metabolomics (ZHU LAB)</p>
      <p>Interdisciplinary Research Center on Biology and Chemistry (IRCBC)</p>
      <p>Shanghai Institute of Organic Chemistry (SIOC)</p>
      <p>Chinese Academy of Sciences (CAS)</p>
      <p>26 Qiuyue Road, Pudong, Shanghai, China 201210</p>
      <p>Website: <a href="http://www.zhulab.cn">www.zhulab.cn</a></p>
    </div>
  </body>
</html>
"""


def dic2string(para_dic):
    """
    convert a python dict to a string as the input of R package parameters
    :param para_dic: a dictionary
    :return: a string
    """
    para_str = ''
    non_string_paras = ['correct', 'p.cutoff', 'threads']
    for key in para_dic.keys():
        if key not in ['control.group', 'case.group']:
            if key in non_string_paras:
                para_str += str(key) + '=' + str(para_dic[key]) + ', '
            # elif key == 'p.cutoff':
            #     para_str += str(key) + '='
            else:
                para_str += str(key) + '="' + str(para_dic[key]) + '", '
    groups = (para_dic['control.group'], para_dic['case.group'])
    para_str += 'group=c("{}", "{}")'.format(groups[0], groups[1])
    return para_str


def get_url(root_dir, pol, result_fn="results.tar.gz",
            log_fn="run.log.txt"):
    """
    Get the urls of result files for website download (web url) and
    the path of log file (local path)
    :param root_dir: MEDIA_ROOT/user_id/project_id (local path), eg: "/tmp/metdna-upload/ffc3-49f4-4559-9979/dd4f85b9"
    :param pol: Positive / Negative / Both
    :param result_fn: result file name, default name is 'results.tar.gz'
    :param log_fn: log file name, default name is 'run.log.txt'
    :return: result url and log path
    @example:
    result url, http://metdna.zhulab.cn/media/a984-1a23-4825-a3f3/fece8f9c/results/results_pos.tar.gz
    log path, /mnt/data/metdna-upload/a984-1a23-4825-a3f3/fece8f9c/POS and NEG/run.log.txt
    """
    pol = pol.lower()
    web_root = root_dir.replace(MEDIA_ROOT, URL)
    result_url_pos = os.path.join(web_root, 'POS', result_fn.replace('results', 'results_pos'))
    result_url_neg = os.path.join(web_root, 'NEG', result_fn.replace('results', 'results_neg'))
    result_url_both = os.path.join(web_root, 'POS and NEG', result_fn.replace('results', 'results_both'))

    log_path_pos = os.path.join(root_dir, 'POS', log_fn)
    log_path_neg = os.path.join(root_dir, 'NEG', log_fn)
    log_path_both = os.path.join(root_dir, 'POS and NEG', log_fn)
    if pol == 'positive':
        return {'result_url': result_url_pos, 'log_path': log_path_pos, 'dl_file_name': result_url_pos.split('/')[-1]}
    elif pol == 'negative':
        return {'result_url': result_url_neg, 'log_path': log_path_neg, 'dl_file_name': result_url_neg.split('/')[-1]}
    else:  # both mode has three parts of result, finally merged them all together
        return {'result_url': ', '.join([result_url_pos, result_url_neg, result_url_both]),
                'log_path': log_path_both, 'dl_file_name': result_url_both.split('/')[-1]}


def check_log_file(url, project_name):
    """
    check log file to confirm the result of MetDNA,
    compress all result files if analysis done successfully
    Also check if the result file exist and delete intermediate data before compressing result
    :param url: the return of function get_url()
    :param project_name: current project name
    contains log_path: a string, '.../POS/run.log.txt' for positive mode; './POS and NEG/run.log.txt' for both mode
                         web url: eg: 'http://localhost:8024/media/ffc3-49f4-4559-9979/7bb4839e/POS/results.tar.gz'
    :return: result_status (done: all is well / error: all is error / pos_error / neg_error)
             exist_status (exist: all is well / error: one or both are not exist)
    """
    # pos: 1, neg: 2
    result_status = 'error'
    exist_status = 'exist'
    log_path = url['log_path']
    result_url = url['result_url'].split(', ')  # a list
    dl_file_name = url['dl_file_name']  # the file name of download result
    print('# result url: ', result_url)
    tar_command_with_log = 'tar -czf {} */ MetDNA.parameters.csv run.log.txt'
    tar_command_no_log = 'tar -czf {} */'
    print('#--------------, in check log file')
    print(tar_command_with_log)
    raw_dir = os.getcwd()
    # deal with url
    result_p_n = []  # local path and file name  [(), ()]
    for _url in result_url:
        _url = _url.replace(URL, MEDIA_ROOT).split('/')
        fn = _url.pop()
        result_p_n.append(('/'.join(_url), fn))

    if os.path.exists(log_path):
        _ = result_p_n
        with open(log_path) as f:
            log = f.readlines()
            log = [i.strip() for i in log]
        if '##MetDNA is done.' in log:
            result_status = 'done'
            # if len(_) == 1:  # positive or negative mode
            #     _path = _[0]
            #     cd_des = _path[0]  # cd destination dir
            #     tar_command = tar_command_with_log.format(_path[1])
            #     print('# cd command: ', cd_des)
            #     print('# tar command: ', tar_command)
            #     os.chdir(cd_des)
            #     os.system(tar_command)
            #     print('current dir: ', os.getcwd())
            #     if os.path.exists(os.path.join(_path[0], _path[1])):
            #         shutil.move(_path[1], '../results/{}'.format(_path[1]))
            #     else:
            #         exist_status = 'error'
            # elif len(_) == 3:  # both mode
            for path in _:
                print('###path is', path, '\n')
                cd_des = path[0]  # a full path of result location
                # single polarity mode or the result folder of both mode that has log file
                if len(_) == 1 or ('POS and NEG' in cd_des):
                    tar_cmd = tar_command_with_log.format(path[1])
                else:
                    tar_cmd = tar_command_no_log.format(path[1])
                print(cd_des, tar_cmd)
                os.chdir(cd_des)
                os.system(tar_cmd)
                if os.path.exists(os.path.join(path[0], path[1])):
                    des_path = os.path.join('..', 'results', path[1])
                    if os.path.exists(des_path):
                        os.remove(des_path)
                    shutil.move(path[1], '../results')
                else:
                    exist_status = 'error'

            # add some code here to merge all results into one folder!
            # and delete intermediate data
            os.chdir('../results')
            # now there are three .tar.gz files (results_pos.tar.gz, results_neg.tar.gz and results_both.tar.gz)
            # for both mode
            # or one .tar.gz file (results_pos.tar.gz or results_neg.tar.gz) for single polarity mode
            all_file_name = [i[1] for i in _]
            for name in all_file_name:
                _dir = re.split('\.|_', name)[1].upper()
                if not os.path.exists(_dir):
                    os.mkdir(_dir)
                extract_cmd = 'tar -xzf {} -C {}'.format(name, _dir)
                os.system(extract_cmd)
                # delete original .tar.gz file
                if os.path.exists(name):
                    os.remove(name)
                # delete intermediate data
                imet_dir1 = os.path.join(_dir, 'MRN_annotation_result', 'intermediate_data')
                imet_dir2 = os.path.join(_dir, 'Pathway_enrichment_analysis_result', 'intermediate_data')
                if os.path.exists(imet_dir1):
                    shutil.rmtree(imet_dir1)
                if os.path.exists(imet_dir2):
                    shutil.rmtree(imet_dir2)
            tar_cmd = tar_command_no_log.format(dl_file_name)
            os.system(tar_cmd)
            # try to rename result file, add project name
            dl_file_name_new = project_name + '_' + dl_file_name
            try:
                os.rename(dl_file_name, dl_file_name_new)
                if os.path.exists(dl_file_name_new):
                    url['dl_file_name'] = dl_file_name_new
            except:
                print('####### mv error ########')
                # print(mv_cmd)
                print(url)

    os.chdir(raw_dir)
    return {'result_status': result_status,
            'exist_status': exist_status,
            'url': url}


# para = 'species="hsa", ce="10", polarity="positive", ' \
#        'neg.path="D://tmp//metdna-upload//22cd-7466-4355-ac40//7db0829a//NEG", sample.info.pos="sample.info.csv", ' \
#        'pos.path="D://tmp//metdna-upload//22cd-7466-4355-ac40//7db0829a//POS", sample.info.neg="sample.info.csv", ' \
#        'uni.test="t", correct=FALSE, ms1.data.neg="", p.cutoff=0.05, column="hilic", ms1.data.pos="data.csv", ' \
#        'group=c("E30", "W15"), threads=10'
# print(get_url(para, file_name='results.zip'))
# print(dic2string(a_dic))
