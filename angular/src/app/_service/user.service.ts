// created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter

import { Injectable } from '@angular/core';
import { Headers, Http, RequestOptions, Response } from "@angular/http";
import { User } from '../_data_model/index';
import { AuthenticationService } from './authentication.service';
import { userApi } from '../_data_model/index';

import 'rxjs/add/operator/toPromise';

@Injectable()
export class UserService {
  private userApi = userApi;

  constructor(
    private http: Http,
    private authService: AuthenticationService,
  ) { }

  getUserByEmail(email: string): Promise<Response> {
    const url = this.userApi
    // console.log('here is in user.service.ts', url);
    let headers = new Headers();
    headers.append('Content-Type', 'application/json');
    let options = new RequestOptions({ headers: headers, params: {'email': email} });
    return this.http.get(url, options)
               .toPromise()
              //  .then(res => res.json().results as User)
               .then(res => res)
              //  .then(res => console.log(res.json().results))
               .catch(this.handleError);
  }
  
  createUser(email: string, name: string, pw: string): Promise<Response> {
    const url = this.userApi;
    let headers = new Headers();
    headers.append('Content-Type', 'application/json');
    let options = new RequestOptions({ headers: headers });
    // console.log(name.length);
    // console.log(name);
    return this.http
            .post(url, JSON.stringify({username: name, 
                                       password: pw, 
                                       email: email,
                                       is_active: true}), options)
            .toPromise()
            // 下面的.then方法对默认返回的数据进行了加工，得到了一个完整的User对象
            .then(res => res)
            .catch(this.handleError);
  }

  userLogin(username: string, pw: string, email: string) {
    return this.authService.login(username, pw, email);
  }

  checkUserToken(token: string) {
    return this.authService.tokenVerify(token);
  }

  userLogout() {
    this.authService.logout();
  }

  getCurrentUser() {
    // {email: '', token: ''}
    return this.authService.getCurrentUser();
  }

  private handleError(error: any): Promise<any> {
    console.error('An error occurred', error); // for demo purposes only
    return Promise.reject(error.message || error);
  }

  // https://stackoverflow.com/a/2117523/2803344
  generateUserId = function() {
    return 'xxxx-xxxx-4xxx-yxxx'.replace(/[xy]/g, function(c) {
      var r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
      return v.toString(16);
    });
  }

}
