// created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter -->

const apiSite: string = 'http://192.168.201.33:6001/';  // 6001 for back-end of formal product in local webserver
// const apiSite: string = 'http://192.168.201.33:6002/';  // 6002 for back-end of test
// const apiSite: string = 'http://localhost:6002/';  // 6002 for back-end of local dev
// const apiSite: string = 'http://metdna.zhulab.cn:6001/';  // the back-end of formal product in cloud webserver

export const frontPort: string = '80';  // 80 for frontend of formal product 
// export const frontPort: string = '8024';  // 8024 for frontend of test

export const tokenAuthApi: string = apiSite + 'api-token-auth/';
export const tokenVerifyApi: string = apiSite + 'api-token-verify/';
export const tokenRefreshApi: string = apiSite + 'api-token-refresh/';
export const userApi: string = apiSite + 'users/';
export const fileApi: string = apiSite + 'files/';
export const projectApi: string = apiSite + 'projects/';
export const paraApi: string = apiSite + 'project-queues/';

