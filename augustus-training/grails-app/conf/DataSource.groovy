dataSource {
}
hibernate {
    cache.use_second_level_cache=true
    cache.use_query_cache=true
    cache.provider_class='com.opensymphony.oscache.hibernate.OSCacheProvider'
}
// environment specific settings
environments {
	development {
		dataSource {
			poold = true
			dbCreate = "update"
			url = "jdbc:mysql://localhost/augtraindev"
                        driverClassName = "com.mysql.jdbc.Driver"
                        username = "grails1"
                        password = "9M6T5kgU"
                        properties {
                                   maxActive = 100
                                   maxIdle = 25
                                   minIdle = 5
                                   initialSize = 5
                                   minEvictableIdleTimeMillis = 60000
                                   timeBetweenEvictionRunsMillis = 60000
                                   maxWait = 10000
                        }                       


		}
	}
	test {
		dataSource {
			poold = true
			dbCreate = "update"
			url = "jdbc:mysql://localhost/augtrainprod"
                        driverClassName = "com.mysql.jdbc.Driver"
                        username = "grails1"
                        password = "9M6T5kgU"
                        properties {
                                   maxActive = 100
                                   maxIdle = 25
                                   minIdle = 5
                                   initialSize = 5
                                   minEvictableIdleTimeMillis = 60000
                                   timeBetweenEvictionRunsMillis = 60000
                                   maxWait = 10000
                        }                       
		}
	}
	production {
		dataSource {
			poold = true
			dbCreate = "update"
			url = "jdbc:mysql://localhost/augtrainprod"
                        driverClassName = "com.mysql.jdbc.Driver"
                        username = "grails1"
                        password = "9M6T5kgU"
                        properties {
                                   maxActive = 100
                                   maxIdle = 25
                                   minIdle = 5
                                   initialSize = 5
                                   minEvictableIdleTimeMillis = 60000
                                   timeBetweenEvictionRunsMillis = 60000
                                   maxWait = 10000
                        }                       
		}
	}
}
