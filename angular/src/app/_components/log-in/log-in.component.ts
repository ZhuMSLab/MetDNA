// created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter

import { Component, OnInit } from '@angular/core';
import {FormControl, Validators} from '@angular/forms';
import { UserService, AuthenticationService } from '../../_service/index';
import { User } from '../../_data_model/index';
import { useAnimation } from '@angular/animations/src/animation_metadata';

const EMAIL_REGEX = /^[a-zA-Z0-9.!#$%&’*+/=?^_`{|}~-]+@[a-zA-Z0-9-]+(?:\.[a-zA-Z0-9-]+)*$/;

@Component({
  selector: 'app-log-in',
  templateUrl: './log-in.component.html',
  styleUrls: ['./log-in.component.css']
})
export class LogInComponent implements OnInit {

  private username: string = '';
  private email: string = '';
  private hideEmailForm: boolean = false;
  public welcomeWords: string = '';
  public logOutInfo: string = '';
  private currentUser: object = {};
  private defaultPW: string = 'fbdd8690';
  private intervalId: any;
  public logInError: string = '';
  private logInErrorMessage: string = 'Can\'t connect to web server, ' +
          'please check your network connection or report this error to us by email.';

  constructor(
    private userService: UserService,
    // private authService: AuthenticationService,
  ) { }

  ngOnInit() {
    // console.log(this.username);
    // let email = localStorage.getItem('email');
    let user = this.userService.getCurrentUser();
    if (user) {
      this.currentUser = user;
      this.email = this.currentUser['email'];
      this.hideEmailForm = true;
    }
    // auto logout after 10 hours of login
    setTimeout(() => {
      this.logOutInfo = 'Signature has expired, please input your email again.';
      this.logout();
    }, 1000*60*60*10);  // expire after 10 hours

    this.intervalId = setInterval(() => {
      this.currentUserStatusMonitor();
    }, 1000 * 60 * 2);  // check user status every 2 minutes
  }

  currentUserStatusMonitor() {
    this.currentUser = this.userService.getCurrentUser();
    // console.log(this.currentUser);
    if (this.currentUser===null) {
      this.logOutInfo = 'Signature has expired, please input your email again.';
      this.logout();
    }
  }

  // https://stackoverflow.com/a/37116635/2803344
  ngOnDestroy() {
    if (this.intervalId) {
      clearInterval(this.intervalId);
    }
  }
  
  emailFormControl = new FormControl('', [
    Validators.required,
    Validators.pattern(EMAIL_REGEX)]);

  onClick() {
    this.username = this.userService.generateUserId();
    this.email = this.emailFormControl.value;
    this.logOutInfo = '';
    
    // console.log(this.username); 
    // console.log(this.email);
    // using getUser before using createUse!!!
    let queryInfo = this.userService.getUserByEmail(this.email);

    queryInfo.then(res => {
      // console.log(typeof(res));
      // console.log(res);
      // console.log(res.json().results);
      let user: User[] = res.json().results;
      if (res.status===200) {
        if (user.length !== 0) {
          console.log('here');
          this.username = user[0].username;
          this.welcomeWords = 'Welcome back!';
          let loginInfo = this.userService.userLogin(this.username, this.defaultPW, this.email);
          console.log(loginInfo);
          loginInfo.then(res => {
            // console.log(res);
            if (res) {
              // get token successfully
              setTimeout(() => {
                this.hideEmailForm = true;
              }, 200)
            }
            else {
              this.logInError = this.logInErrorMessage;
            }
          });
        }
        else {
          // create a new user
          let createInfo = this.userService.createUser(this.email, this.username, this.defaultPW);
          createInfo.then(res => {
            // console.log(res.json()[0]);  // we can get returned user like this
            // console.log(this.email);
            if (res.status === 201) {  // created
              this.welcomeWords = 'Create new account successfully!';
              console.log(this.welcomeWords);
              let loginInfo = this.userService.userLogin(this.username, this.defaultPW, this.email);
              console.log(loginInfo);
              loginInfo.then(res => {
                // console.log(res);
                if (res) {
                  // get token successfully
                  setTimeout(() => {
                    this.hideEmailForm = true;
                  }, 200)
                }
                else {
                  this.logInError = this.logInErrorMessage;
                }
              });
            }
            else {
              this.logInError = this.logInErrorMessage;
            }
          }).catch(res => {
            this.logInError = this.logInErrorMessage;
          });
        }
      }
      else {
        this.logInError = this.logInErrorMessage;
      }
    }).catch(res => {
      this.logInError = this.logInErrorMessage;
    });
  }

  logout() {
    this.userService.userLogout();
    this.emailFormControl.reset();
    this.currentUser = {};
    this.email = '';
    this.hideEmailForm = false;
    this.welcomeWords = '';
    this.logInError = '';
  }
}
