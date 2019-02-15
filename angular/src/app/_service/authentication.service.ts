// created by Xin Xiong<onlybelter@outlook.com>
import { Injectable } from '@angular/core';
import { Http, Headers, Response, RequestOptions } from '@angular/http';
import { Observable } from 'rxjs';
import { tokenAuthApi, tokenRefreshApi, tokenVerifyApi } from '../_data_model/index';
import 'rxjs/add/operator/map';
import 'rxjs/add/operator/timeout';

@Injectable()
export class AuthenticationService {
  public token: string;
  private url_auth = tokenAuthApi;
  private url_verify = tokenVerifyApi;
  private url_refresh = tokenRefreshApi;
  // private options: RequestOptions;
  private currentUser: any;

  constructor(
    private http: Http,
  ) { }

  login(username: string, password: string, email: string): Promise<boolean|number> {
    // console.log('I am in login');
    // console.log(username);
    // console.log(password);
    let headers = new Headers();
    headers.append('Content-Type', 'application/json');
    let options = new RequestOptions({ headers: headers });
    return this.http.post(this.url_auth, 
                          JSON.stringify({ username: username, password: password }), options)
                          .toPromise()
                          .then(res => {
                            // console.log(res);
                            let token = res.json() && res.json().token;
                            if (token) {
                              // store username and jwt token in local storage to keep user logged in between page refreshes
                              localStorage.setItem('currentUser', 
                                                   JSON.stringify({ email: email, token: token }));
                              return true;
                            }
                            else {
                              return false;
                            }
                          })
                          .catch(this.handleError);
  }

  logout(): void {
    //remove user from local storage to log user out
    localStorage.clear();

  }

  getCurrentUser(): Observable<any> {
    this.currentUser = JSON.parse(localStorage.getItem('currentUser'));
    return this.currentUser;
  }

  tokenVerify(token: string): Promise<number> {
    let headers = new Headers();
    headers.append('Content-Type', 'application/json');
    let options = new RequestOptions({ headers: headers });
    return this.http.post(this.url_verify, JSON.stringify({'token': token}), options)
            .toPromise()
            .then((response: Response) => {
              // console.log(response);
              return response.status;
            })
            .catch(this.handleError);
  }

  private handleError(error: any): Promise<number> {
    // console.error('An error occurred', error);
    // console.log(error.status);
    return Observable.of(error.status).toPromise();
  }
}

